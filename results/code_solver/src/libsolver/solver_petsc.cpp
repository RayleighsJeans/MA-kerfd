#include <iostream>
#include <vector>
#include <string>
#include <fstream>

#include "epsilon.h"
#include "use_flags.h"
#include "grid.h"
#include "var.h"
#include "solver.h"
#include "fill_rhs.h"
#include "mpi_wrapper.h"
#include "benchmark.hpp"

#if !NTRL_ONLY
#include <petsc.h>
#include <petscksp.h>
#include <mpi.h>
#endif

#define DEBUG_PETSC 0

namespace{
#if !NTRL_ONLY
  static char help[] = "Solves linear system with KSP.\n\n";
  Vec phi_petsc, rhs_petsc; //approx solution, RHS
  Mat matrix; //linear system matrix
  KSP ksp; //linear solver context
  PC pc; //preconditioner context
  PetscMPIInt    petsc_size; //PETSC_COMM_WORLD communicator size
  DM             da; //distributed array
  DMDALocalInfo  info; //local sizes within dmda
  KSPConvergedReason reason; //reason why a method did/didnt converge
  PetscLogStage  stage; //for logging information
  PetscViewer    viewer; //petsc viewer object
  std::vector<double> dens; //necessary for fill_rhs
#endif
}

//because types of petsc functions dont match...
PetscErrorCode KSPMonitorTrueResidualNormWrapper(KSP ksp, PetscInt n, PetscReal rnorm, void* dummy){
  KSPMonitorTrueResidualNorm(ksp,n,rnorm,(PetscViewerAndFormat*)dummy);
  return 1;
}

void init_solver(std::string filename, int argc, char** argv){
#if !NTRL_ONLY
  //global size parameters
  int NR = global_grid.r_dim+1;
  int NZ = global_grid.z_dim+1;
  int NG = NR*NZ;

  //set dens and eps arrays to necessary sizes
  dens.resize(NG,0);

  //calculate epsilon and read from file
  std::vector<std::vector<double>> eps; //vector for epsilon matrix
  set_solver_boundary(filename, eps); 

  //init Petsc
  PetscInitialize(&argc,&argv,(char*)0,help);

#if DEBUG_PETSC
  //allow verbose information to be printed
  PetscInfoAllow(PETSC_TRUE,"petsc_info.log");

  //set global format for VecView
  PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_NATIVE);
#endif
  
  //PETSC_COMM_WORLD defaults to MPI_COMM_WORLD if MPI is already up and running
  MPI_Comm_size(PETSC_COMM_WORLD,&petsc_size);
  
  //create context for distributed rectangular 2d arrays (degr. of freedom=1, stencil_width=1=star-shaped 5pt-stencil, PETSC_DECIDE for local dimensions
  //calculate number of processors and grid points in each dimension
  {
    int r=1;
    int z=1;
    int mpi_size = get_size();
    for(auto i=0;r*z<mpi_size;++i){
      if(NZ/z < NR/r) z++;
      else r++;
    }
    DMDACreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_GHOSTED,DM_BOUNDARY_GHOSTED,DMDA_STENCIL_STAR,NR,NZ,r,z,1,1,NULL,NULL,&da);
  }

  //retrieve local size informations
  DMDAGetLocalInfo(da,&info);

  //create parallel DMDA matrix
  DMCreateMatrix(da,&matrix);

  //matrix assembly, set a petsc log stage string and then assemble matrix for local dimensions 
  //Elements belonging to other processe will be sent during MatAssembly
#if DEBUG_PETSC
  std::cout<<"PETSC local dimensions on rank "<<get_rank()<<": rmin="<<info.xs<<" rmax="<<info.xs+info.xm<<" zim="<<info.ys<<" zmax="<<info.ys+info.ym<<std::endl;
  PetscLogStageRegister("Matrix Assembly",&stage);
  PetscLogStagePush(stage);
#endif
  for(auto r=info.xs;r<info.xs+info.xm;++r){
    for(auto z=info.ys;z<info.ys+info.ym;++z){
      int ind = r*NZ+z;		 
      MatStencil row[1];
      row[0].i = r; row[0].j = z;
        
      if ( eps[r][z] == -100. ) // Fixed potential on metal
      {
	MatStencil col[1]; 
	PetscScalar val[1];
	//bingo
	col[0].i = r; col[0].j = z;
	val[0] = 1;

        MatSetValuesStencil(matrix,1,row,1,col,val,INSERT_VALUES);
      } 
      else if ( eps[r][z] == -200. ) // Fixed gradient rhs
      {
        MatStencil col[2];
        PetscScalar val[2];
	//west
	col[0].i = r; col[0].j = z-1;
	val[0]=-1;
	//bingo
	col[1].i = r; col[1].j = z;
	val[1]=1;

        MatSetValuesStencil(matrix,1,row,2,col,val,INSERT_VALUES);
      } 
      else if ( eps[r][z] == -300. ) // Fixed gradient lhs 
      {
        MatStencil col[2];
        PetscScalar val[2];
	//bingo
	col[0].i = r; col[0].j = z;
	val[0]=1;
	//east
	col[1].i = r; col[1].j = z+1;
	val[1]=-1;

        MatSetValuesStencil(matrix,1,row,2,col,val,INSERT_VALUES);
      } 
      else if ( eps[r][z] == -400. ) // Fixed gradient upper side  
      {
        MatStencil col[2];
        PetscScalar val[2];
	//south
	col[0].i = r-1; col[0].j = z;
	val[0]=-1;
	//bingo
	col[1].i = r; col[1].j = z;
	val[1]=1;

        MatSetValuesStencil(matrix,1,row,2,col,val,INSERT_VALUES);
      } 
      else if ( r == 0)          // Axis 
      {
        //calc coefficients
        double e_iP1_jM1 = eps[r+1][z-1];
        if ( e_iP1_jM1 <= 0.)   e_iP1_jM1 = eps[r][z];
        
        double e_iP1_jP1 = eps[r+1][z+1];
        if ( e_iP1_jP1 <= 0.)   e_iP1_jP1 =  eps[r][z];  
        
        double west = e_iP1_jM1;                                                 // a,b,c,d,e    *= dz*dz
        double east = e_iP1_jP1;
        double north = 2.*(e_iP1_jM1 + e_iP1_jP1);
        double bingo = -3.*(e_iP1_jM1 + e_iP1_jP1);    

        MatStencil col[4];
        PetscScalar val[4];

	//west
	col[0].i = r; col[0].j = z-1;
        val[0]=west;

	//bingo
	col[1].i = r; col[1].j = z;
        val[1]=bingo;

	//east
	col[2].i = r; col[2].j = z+1;
        val[2]=east;

	//north
	col[3].i = r+1; col[3].j = z;
        val[3]=north;

        MatSetValuesStencil(matrix,1,row,4,col,val,INSERT_VALUES);
      }
      else              // the rest: channell, dielectric, boundaries, etc           
      {
        double e_iP1_jM1 =  eps[r+1][z-1];
        if ( e_iP1_jM1 <= 0.)   e_iP1_jM1 = eps[r][z];
        if ( e_iP1_jM1 ==0 )
        { 
          if (  eps[r+1][z]  > 0 ) e_iP1_jM1 = eps[r+1][z];
          else                           e_iP1_jM1 = eps[r][z-1];   
        }
                    
        double e_iP1_jP1 = eps[r+1][z+1];
        if ( e_iP1_jP1 <= 0.)   e_iP1_jP1 = eps[r][z];  
        if ( e_iP1_jP1 ==0 ) 
        {
          if (  eps[r+1][z]  > 0 )   e_iP1_jP1 = eps[r+1][z];
          else      e_iP1_jP1 = eps[r][z+1];   
        }
       
        double e_iM1_jM1 =  eps[r-1][z-1];
        if ( e_iM1_jM1 <= 0.)   e_iM1_jM1 = eps[r][z];
        if (e_iM1_jM1 ==0 ) 
        {
          if ( eps[r-1][z] > 0 ) e_iM1_jM1 = eps[r-1][z];
          else                         e_iM1_jM1 = eps[r][z-1];  
        }      
                        
        double e_iM1_jP1 = eps[r-1][z+1];
        if ( e_iM1_jP1 <= 0.)   e_iM1_jP1 = eps[r][z];
        if (e_iM1_jP1 ==0 )
        { 
          if (  eps[r-1][z] > 0 )  e_iM1_jP1 =  eps[r-1][z];
          else                           e_iM1_jP1 = eps[r][z+1];   
        }
              
        double north = 0.25*(e_iP1_jM1 + e_iP1_jP1)/r + 0.5*(e_iP1_jM1 + e_iP1_jP1);
        double south = 0.5*(e_iM1_jM1 + e_iM1_jP1) -  0.25*(e_iM1_jM1 + e_iM1_jP1)/r;
        double bingo = 0.25*(e_iM1_jM1 + e_iM1_jP1 - e_iP1_jM1 - e_iP1_jP1)/r - (e_iP1_jM1 + e_iP1_jP1 + e_iM1_jM1 + e_iM1_jP1 );
        double west = 0.5*(e_iP1_jM1 + e_iM1_jM1 );
        double east = 0.5*(e_iP1_jP1 + e_iM1_jP1 ); 

        MatStencil col[5];
        PetscScalar val[5];

	//south
	col[0].i = r-1; col[0].j = z;
	val[0]=south;

	//west
	col[1].i = r; col[1].j = z-1;
        val[1]=west;

	//bingo
	col[2].i = r; col[2].j = z;
        val[2]=bingo;

	//east
	col[3].i = r; col[3].j = z+1;
        val[3]=east;

	//north
	col[4].i = r+1; col[4].j = z;
        val[4]=north;

        MatSetValuesStencil(matrix,1,row,5,col,val,INSERT_VALUES);
      }
    }
  }

  MatAssemblyBegin(matrix,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(matrix,MAT_FINAL_ASSEMBLY);
#if DEBUG_PETSC
  PetscLogStagePop();
#endif

  //create global vectors
  DMCreateGlobalVector(da,&phi_petsc);
  VecDuplicate(phi_petsc,&rhs_petsc);

  //create linear solver context and set system and preconditioner matrix
  KSPCreate(PETSC_COMM_WORLD,&ksp);
  KSPSetOperators(ksp,matrix,matrix);

  /***********Configure Options for the solver***********/
  //set the solver type here (KSPGMRES = Generalized Minimal RESidual)
  //best options are (in that order) KSPFGMRES, KSPGCR, KSPGMRES, KSPDGMRES
  KSPSetType(ksp,KSPFGMRES);

  //set preconditioning side (FGMRES only supports PC_RIGHT)
  KSPSetPCSide(ksp,PC_RIGHT);

  //tell ksp that non-zero initial guess is to be used
  KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);

  //configure pre-conditioner
  //extract the pre-conditioner
  KSPGetPC(ksp,&pc);
  //set the preconditioner 
  PCSetType(pc,PCHYPRE);

  /**********************************************************************************************************************/
  {
    //configure hypre boomeramg options here
    //threshold for strong connections (default=0.25)
    PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_strong_threshold","0.25");

    //interpolation truncation factor (default=0=no truncation)
    PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_truncfactor","0.01");

    //set cycle type ("W" does not scale well for parallel applications. Can be {"V","W"} (default=V)
    PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_cycle_type","V");

    //set max number of AMG levels, at least 2 (default=25)
    PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_max_levels","25");
    
    //set max number of iterations per hypre call (default=1)
    PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_max_iter","1");
    
    //set convergence tolerance for hypre call (default=0=fixed number of iterations)
    PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_tol","0");
    
    //maximum number of row elements used in interpolation operator (default=0=unlimited)
    PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_P_max","0");

    //maximum row sum, must be between 0 and 1 (default=0.9)
    PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_max_row_sum","0.9");

    //use nodal coarsening from 1-6. (default=0=no nodal coarsening) (good values=0,6,3)
    PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_nodal_coarsen","0");

    //use vector interpolation variant from 1-3. (default=0=no interp variant)
    PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_vec_interp_variant","0");
    
    //set coarsening type. Can be any of (default=Falgout):
    //{"CLJP","Ruge-Stueben","modifiedRuge-Stueben","Falgout","PMIS"(not working),"HMIS"}
    PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_coarsen_type","HMIS");
    //number of levels/paths for aggressive coarsening (default=0/default=1)
    PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_agg_nl","0");
    PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_agg_num_paths","1");

    //set measure type. Can be {"local","global"} (default=local)
    PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_measure_type","local");

    //more complex smoother type can be any of (default="Schwarz-smoothers"):
    //{"Schwarz-smoothers","Pilut","ParaSails","Euclid"}
    //comment the line if no complex smoother is desired and/or set number of levels to 0
    //PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_smooth_type","Schwarz-smoothers");
    //set number of levels where more complex smoother is used(default=0)
    PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_smooth_num_levels","0");

    //relaxation type for going up/down/coarsest level can be any of (default=symmetric-SOR/Jacobi):
    //{"Jacobi","sequential-Gauss-Seidel","seqboundary-Gauss-Seidel",
    //"SOR/Jacobi","backward-SOR/Jacobi","symmetric-SOR/Jacobi", "l1scaled-SOR/Jacobi",
    //"Gaussian-elimination","CG","Chebyshev","FCF-Jacobi","l1scaled-Jacobi"}
    //PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_relax_type_all","symmetric-SOR/Jacobi");
    PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_relax_type_down","symmetric-SOR/Jacobi");
    PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_relax_type_up","symmetric-SOR/Jacobi");
    PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_relax_type_coarse","Gaussian-elimination");

    //relaxation weight on all levels (default=1=weight everywhere, positive number=weight, 0=hypre estimates, -k=determined by k CG steps)
    PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_relax_weight_all","1");

    //do not use CF-Relaxation (default=0, so CF is used) (1 (with otherwise default params) is 10% faster for HEMP)
    PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_no_CF","1");

    //number of grid sweeps for all levels, or up, down and coarsest (default=1 for all levels)
    //PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_grid_sweeps_all","1");
    PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_grid_sweeps_down","1");
    PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_grid_sweeps_up","2");
    PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_grid_sweeps_coarse","1");

    //interpolation type can be any of (default=classical)
    //{"classical","direct","multipass","multipass-wts","ext+i","ext+i-cc","standard","standard-wts","FF","FF1"}
    PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_interp_type","ext+i-cc");

    //print statistics (default=0 for production runs)
#if DEBUG_PETSC
    PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_print_statistics","1");

    //print debug (default=0 for production runs)
    PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_print_debug","1");
#endif
  }
  /**********************************************************************************************************************/
  
  
  //set parameters to determine wether convergenc ewas achieved
  //PETSC_DEFAULT=standard parameters
  KSPSetTolerances(ksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);

#if DEBUG_PETSC
  KSPMonitorSet(ksp,KSPMonitorTrueResidualNormWrapper,nullptr,nullptr);
#endif

  //set up ksp solver from options set above
  KSPSetFromOptions(ksp);
#endif
}

void destroy_solver(){
#if !NTRL_ONLY
  //Free work space.  All PETSc objects should be destroyed when they
  //are no longer needed.
  dens.resize(0);
  MatDestroy(&matrix);
  VecDestroy(&phi_petsc);
  VecDestroy(&rhs_petsc);

  //destroy contexts and solver
  DMDestroy(&da);
  KSPDestroy(&ksp);

  //Always call PetscFinalize() before exiting a program.  This routine
  //Finalizes the PETSc libraries as well as MPI (if MPI_Init was invoked by petsc)
  //Provides summary and diagnostic information if certain runtime options are chosen (e.g., -log_summary).
  PetscFinalize();
#endif
}

//fill rhs and put it into petsc vector for solution
void calculate_rhs(){
#if !NTRL_ONLY
  int NR,NZ; //global grid size
  PetscScalar** arr; //temp array pointing to rhs values

  //retrieve grid dimensions
  //TODO do this globally in this file
  DMDAGetInfo(da,0,&NR,&NZ,0,0,0,0,0,0,0,0,0,0);

  //get rhs from charge densities
  fill_rhs(dens_el, dens_ion, dens_nion, dens_ion2p, surf_densR, surf_densZ,dens.data(), qe, NR, NZ,nrMAXmodel, nzMINmodel, nzMAXmodel);

  //get temporary 2D array associated with rhs_petsc to store rhs values in
  DMDAVecGetArray(da,rhs_petsc,&arr);

  //petsc uses reverse indexing... 
  for(auto z=info.ys;z<info.ys+info.ym;++z){
    for(auto r=info.xs;r<info.xs+info.xm;++r){
      int ind = r*NZ+z;
      arr[z][r] = dens[ind];
    }
  }

  //return temporary array to petsc voodoo
  DMDAVecRestoreArray(da,rhs_petsc,&arr);

  //communicate values to all other processes
  VecAssemblyBegin(rhs_petsc);
  VecAssemblyEnd(rhs_petsc);
#endif
}

void calculate_potential( double* field_time, int r1, int r2, int z1, int z2 ){
#if !NTRL_ONLY
  auto benchstart = Benchmark::start();
  int NR,NZ; //global array size
  int its = 0; //number of iterations
  PetscScalar** arr; //array to store rhs values in

  //calculate rhs of the system
  calculate_rhs();
  
  //solve linear system 
  KSPSolve(ksp,rhs_petsc,phi_petsc);

#if DEBUG_PETSC
  //retrieve convergence info
  KSPGetConvergedReason(ksp,&reason);

  //check wether solver gave diverging results for debugging purposes
  if (reason==KSP_CONVERGED_RTOL){
    PetscPrintf(PETSC_COMM_WORLD,"\nConvergence reason return value=KSP_CONVERGED_RTOL.\n");
  }else if(reason==KSP_CONVERGED_ATOL){
    PetscPrintf(PETSC_COMM_WORLD,"\nConvergence reason return value=KSP_CONVERGED_ATOL.\n");
  }else if(reason==KSP_CONVERGED_ITS){
    PetscPrintf(PETSC_COMM_WORLD,"\nConvergence reason return value=KSP_CONVERGED_ITS.\n");
  }else if(reason==KSP_CONVERGED_CG_NEG_CURVE){
    PetscPrintf(PETSC_COMM_WORLD,"\nConvergence reason return value=KSP_CONVERGED_CG_NEG_CURVE.\n");
  }else if(reason==KSP_CONVERGED_CG_CONSTRAINED){
    PetscPrintf(PETSC_COMM_WORLD,"\nConvergence reason return value=KSP_CONVERGED_CG_CONSTRAINED.\n");
  }else if(reason==KSP_CONVERGED_STEP_LENGTH){
    PetscPrintf(PETSC_COMM_WORLD,"\nConvergence reason return value=KSP_CONVERGED_STEP_LENGTH.\n");
  }else if(reason==KSP_CONVERGED_ITERATING){
    PetscPrintf(PETSC_COMM_WORLD,"\nConvergence reason return value=KSP_CONVERGED_ITERATING.\n");
  }else if(reason==KSP_DIVERGED_ITS){
    PetscPrintf(PETSC_COMM_WORLD,"\nDivergence reason return value=KSP_DIVERGED_ITS.\n");
  }else if(reason==KSP_DIVERGED_DTOL){
    PetscPrintf(PETSC_COMM_WORLD,"\nDivergence reason return value=KSP_DIVERGED_DTOL.\n");
  }else if(reason==KSP_DIVERGED_NANORINF){
    PetscPrintf(PETSC_COMM_WORLD,"\nDivergence reason return value=KSP_DIVERGED_NANORINF.\n");
  }else if(reason==KSP_DIVERGED_BREAKDOWN){
    PetscPrintf(PETSC_COMM_WORLD,"\nDivergence reason return value=KSP_DIVERGED_BREAKDOWN.\n");
  }else if(reason==KSP_DIVERGED_BREAKDOWN_BICG){
    PetscPrintf(PETSC_COMM_WORLD,"\nDivergence reason return value=KSP_DIVERGED_BREAKDOWN_BICG.\n");
  }else if(reason==KSP_DIVERGED_INDEFINITE_PC){
    PetscPrintf(PETSC_COMM_WORLD,"\nDivergence reason return value=KSP_DIVERGED_INDEFINITE_PC.\n");
  }else{
    PetscPrintf(PETSC_COMM_WORLD,"\nConvergence reason return value unknown=%d.\n",(int) reason);
  }

  //View solver info; we could instead use the option -ksp_view to
  //print this info to the screen at the conclusion of KSPSolve().
  KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);

  /*
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Check solution and clean up
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  
  //Check the error
  KSPGetIterationNumber(ksp,&its);
  PetscPrintf(PETSC_COMM_WORLD,"PETSc_Iterations %D\n",its);
#endif


  //copy solution to code array
  //TODO do this globally before
  DMDAGetInfo(da,0,&NR,&NZ,0,0,0,0,0,0,0,0,0,0);
  DMDAVecGetArray(da,phi_petsc,&arr);
  //zero phi array before this
  //TODO not really efficient
  std::fill(&phi[0],&phi[NG],0);
  for(auto z=info.ys;z<info.ys+info.ym;++z){
    for(auto r=info.xs;r<info.xs+info.xm;++r){
      int ind = r*NZ+z;
      phi[ind] = arr[z][r];
    }
  }
  DMDAVecRestoreArray(da,phi_petsc,&arr);

  //otherwise petsc might apparently cause a deadlock here
  barrier();

  allreduce(&phi[0],NR*NZ);

  long bench_time = Benchmark::stop(benchstart);
  reduce_max(&bench_time);
  *field_time += bench_time;
#endif
}
