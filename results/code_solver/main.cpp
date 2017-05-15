/***********************************************************
2d3v Axially Symmetric PIC-MCC ( ASPIC )

24.01.10

************************************************************/

// this is needed to tell the compiler 
// to define the variables marked with extern in this translation
// unit. ALL OTHER FILES MUST NOT SET THIS FLAG
#define DEFINE_VARIABLES 

#include <limits.h>
#include <sys/stat.h>
#include <string>
#include <iostream>
#include <fstream>

#include "pic.h"
#include "h5save.h"
#include "react.h"
#include "arrays.h"
#include "var_dim.h"

#include "macros.h"

#include "mpi_wrapper.h"

#include "bfield.h"
#include "init.h"
#include "boundary_check.h"
#include "init_prtcls.h"
#include "fill_rhs.h"
#include "collisions.h"
#include "push.h"

#include "efield.h"
#include "read.h"
#include "output.h"
#include "fem_solver.h"
#include "moms.h"
#include "memory.h"
#include "init_diag.h"
#include "save.h"
#include "e_ion.h"
#include "engy.h"
#include "aver.h"
#include "epsilon.h"
#include "emission.h"
#include "domain_decomposition.h"
#include "diag_flags.h"
#include "restart.h"
#include "backup_and_restore.h"
#include "injection.h"
#include "solver.h"
#include "use_flags.h"
#include "debug_flags.h"
#include "domain_init.h"
#include "particle_tracing.h"
#include "benchmark.hpp"
#include "grid.h"
#include "signal_handling.h"
#include "program_options.h"

#define DEBUG_LEVEL DEBUG_ERROR  //DEBUG_INFO_3
#include "debug_printing.h"

#include "density.h"
#include <unistd.h>
#include <math.h>
#include <getopt.h>
#include <limits>

#define exists(file) (stat(file, &Statbuf) < 0 ? 0 : Statbuf.st_mode)

// check if therer are some inf or nan's -> treminate the run
// TODO: check if it limits performance
#include <fenv.h>

struct stat Statbuf;
//move to var.h or another file?
int switch_backup_file=0;

/// @brief deletes all particles that are dont fit into this domain 
void filter_particles_by_grid_size(GridLayer& layer) 
{
  layer.for_all_cells([&layer](Cell& cell){
      std::vector<Particle>& vp = cell.particles;
      for(unsigned int i=0;i<vp.size();++i){
	Particle& pt=cell.particles[i];
	if ( boundary_grid( pt, layer.r_dim, layer.z_dim ) ) cell.remove(i);
      }
    }
  );
}

__attribute__((noreturn))
void terminate_execution( 
			  double field_time, 
			  double pushE_time, 
			  double coll_time, 
			  double injN_time,
			  double pushN_time,
			  double colN_time,
			  double inj_time,
			  double order_time,
			  double orderN_time,
			  Bench_Start begin_main_time,
			  double engcheck,
			  Vec3d momcheck
){
  int mpi_rank = get_rank();
  printf("save at the END: nstep= %i \n",nstep);

  switch_backup_file = !switch_backup_file;
#if GEOM_TYPE_HEMP
#if !NTRL_ONLY
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      //TEMPORARY HOTFIX FOR SURFACE DENSITY SAVE
      //DELETE ONCE CELL SAVE STRUCTURE IS INTRODUCED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      for(int i=0; i<global_grid.z_dim;i++){
	surf_densZ[i]=0;
      }
      for(int specie=0; specie<MAX_SPECIES;specie++){
	GridLayer& layer=global_grid.layers[specie];
	for(int j=0; j<global_grid.z_dim;j++){
	  Cell& cell=layer.get_cell(nROUTthr, j);
	  surf_densZ[j]+=cell.surf_dens[1];
	  surf_densZ[j+1]+=cell.surf_dens[2];
	}
      }
      reduce(surf_densZ, global_grid.z_dim+1);
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif 
#endif 
  std::string filename = "save/particle_backup" + std::to_string(switch_backup_file) + ".h5";
  h5_backup_particle_data(filename.c_str());

  printf("save particles finished\n");

  /// free memory 
  //printf("now free memeory:\n");
  free_arrays();
  //printf("now destroy solver stuff:\n");
  destroy_solver();

  if ( mpi_rank == 0 ){
    printf(
	"Solver: %11.3e   PushE: %11.3e   Coll: %11.3e   Inj: "
	"%11.3e    PushN: %11.3e    CollN: %11.3e   InjN: %11.3e [ms]\n",
	field_time, pushE_time, coll_time, inj_time, pushN_time, colN_time, injN_time);
    printf(
	" EngCheck = %e  MomCheck.r = %e  "
	" MomCheck.t = %e  MomCheck.z = %e\n \n",
	engcheck, momcheck.r, momcheck.t, momcheck.z);
    double main_time = Benchmark::stop(begin_main_time);
    printf("Main_Time  %e ms \n", main_time );   
  }

  finalize_mpi();
  exit(EXIT_SUCCESS);
}

int main(int argc, char** argv) {

#if BUILD_DEBUG
    // check if therer are some inf or nan's -> terminate the run
    feenableexcept(FE_INVALID | FE_OVERFLOW);
#endif


    ProgramOptions program_options;
    int mpi_rank = 0;
    printf("init mpi\n");
    initialize_mpi(&argc,&argv,mpi_rank);
    print0("init signal handling\n");
    init_signal_handling();
    print0("init program options\n");
    init_program_options(argc, argv, program_options);
    print0("init parameters\n");
    init_parameters();
    print0("init random numbers\n");
    seed = 1234 + mpi_rank;
    init_rand(&seed, 1, 52, 0, 0);
    print0("init arrays\n");
    init_arrays();
    print0("init bfield\n");
    init_bfield();
    print0("init reactions\n");
    init_reactions();
    print0("init global grid\n");
    init_global_grid();
    print0("init diag\n");
    init_diag();
    print0("init solver\n");
    init_solver(geometry_file,argc,argv);
    print0("flags used in this project: \n");
    if ( get_rank() == 0 ){
      // standard output stream (terminal/slurm)
      PRINT_USE_FLAGS(std::cout);
      PRINT_DIAG_FLAGS(std::cout);
      PRINT_DEBUG_FLAGS(std::cout);
      // output of flags into file
      std::ofstream flags_output( "FLAGS_ALL.dat", std::ios::app );
      PRINT_USE_FLAGS( flags_output ); 
      PRINT_DIAG_FLAGS( flags_output );
      PRINT_DEBUG_FLAGS( flags_output );
    }
 
    // let all processors read all particle data and afterwards filter everything
    if (program_options.restart) {
      h5_restore_particle_data("particle_backup.h5");

#if GEOM_TYPE_HEMP
#if !NTRL_ONLY
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      //TEMPORARY HOTFIX FOR SURFACE DENSITY READ FROM OLD RUNS TO MAP ONTO CELLS
      //DELETE ONCE CELL SAVE STRUCTURE IS INTRODUCED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(get_rank()==0){
	std::cout<<"Adding loaded surface density charge in cells at r="<<nROUTthr<<std::endl;
	for(int j=0; j<global_grid.z_dim;j++){
	  Cell& cell=global_grid.layers[ELECTRONS].get_cell(nROUTthr, j);
	  cell.surf_dens[1]=+surf_densZ[j];
	}
      }
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif 
#endif 

#if GEOM_TYPE_EADS
#if !NTRL_ONLY
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      //TEMPORARY HOTFIX FOR SURFACE DENSITY READ FROM OLD RUNS TO MAP ONTO CELLS
      //DELETE ONCE CELL SAVE STRUCTURE IS INTRODUCED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(get_rank()==0){
	std::cout<<"Adding loaded surface density charge in cells at r="<<nROUTthr<<std::endl;
	for(int i=0; i<global_grid.r_dim;i++){
	  Cell& cell=global_grid.layers[ELECTRONS].get_cell(i, nzMINmodel);
	  cell.surf_dens[1]=+surf_densR[i];
	}
	for(int j=0; j<global_grid.z_dim;j++){
	  Cell& cell=global_grid.layers[ELECTRONS].get_cell(nrMAXmodel, j);
	  cell.surf_dens[1]=+surf_densZ[j];
	}
      }
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif 
#endif 

      // TODO at the point we have c++ put this into the constructor 
#if USE_PARTICLE_WEIGHTING
      // TODO this wont work anymore
      for (int i = 0; i < ne; ++i){ electrons[i].w = 1; }
      for (int i = 0; i < ni; ++i){ ions[i].w = 1; }
      for (int i = 0; i < nn; ++i){ neutrals[i].w = coll_fac_ntrl_ntrl; }
#if USE_FEEDGAS_OXYGEN
      for (int i = 0; i < nin; ++i){ nions[i].w = 1; }
#endif
#endif
      if (program_options.diag==true && mpi_rank == 0) {
	restore_diagnostic("save/diagnostic.dat");
      }
      nav_start = nstep;
      minstep = nstep + 1;
    }else{
      if(mpi_rank==0){
	//if not restart, inject particles here only on rank 0
	injection_init();
      }
    }

    // delete all particles which are out of the current domain
    print0("filter particle by grid size: nr=%i  nz=%i\n",nr,nz);
    for (int i = 0; i < MAX_SPECIES; ++i){
      filter_particles_by_grid_size( global_grid.layers[i] );
    }

    global_grid.for_all_layers( [](GridLayer& layer){ 
      print0("%s=%zd ",layer.name.c_str(), layer.size());
    } );
    //printf("(rank=%i)\n",mpi_rank);

    //order all particles before init_domain_decomposition here
    //all initial particle injections need to be done before this point
    //global_grid.order_all();

    init_domain_decomposition();

    // TODO this initialization stuff can be done in init_solver 
    int r1 = 0 ,r2 = global_grid.mesh_r_dim;
    int z1 = 0, z2 = global_grid.mesh_z_dim;

    //get domain bounds for each core
    get_domain_bounds(&r1,&r2,&z1,&z2);

    decompose_domain_all( global_grid );

    global_grid.for_all_layers( [](GridLayer& layer){ 
      if ( layer.specie == NEUTRALS ) return;
      density( layer ); 
    });
   
    //printf("fill rhs (rank=%i)\n",mpi_rank);
    fill_rhs(phi, qe);

    if (!program_options.restart) {
	ionstep = dt_ion;
	ntrlstep = dt_ntrl;
	minstep = 1;
      } 

      int steps_offset = atoi(program_options.str_steps);
      if ( steps_offset != -1 ){
	maxstep = minstep + steps_offset;
      }else{
	maxstep = std::numeric_limits<int>::max();
      }

    global_grid.print_size_all();
    print0("nstep=%d minstep=%d maxstep=%d ionstep=%d ntrlstep=%d rank=%d\n",nstep,minstep,maxstep,ionstep,ntrlstep,mpi_rank);
    print0("calculating from %d to %d rank=%d\n",minstep,maxstep-1,mpi_rank);
    print0("first output at %d, then every %d rank=%d\n",nav_start+nav_time,nav_time,mpi_rank);

    // time measurement of the benchmark
    double  field_time = 0., pushE_time = 0., coll_time = 0., 
	    injN_time = 0., pushN_time = 0., colN_time = 0., 
	    orderN_time = 0, order_time = 0., inj_time = 0.; 
    
    if ( get_rank() == 0 ) {
    	std::ofstream output;
	output.open( "benchmark.csv");
    	output << "#nstep,field_time,pushE_time,coll_time,injN_time,"; 
    	output << "pushN_time,colN_time,Inject_time" << std::endl;
    }
    
    double engcheck = 0.;
    Vec3d momcheck;
    momcheck.r = momcheck.t = momcheck.z = 0.0;

#if !NTRL_ONLY
    calculate_potential( &field_time, r1,r2,z1,z2 ); 
    //TODO @ all: this creates slightly increased E_ion for first subcycling after startup!
    electric_field(phi, E_grid, E_ion);
#endif

    if ( mpi_rank != 0 ){
      // reset all diagnotics
      global_grid.for_all_cells( [](Cell& cell){ cell.clear_diagnostics(); } );
    }

#if NTRL_ONLY
    print0("------------------------------------------------------------\n");
    print0("This is a neutral only run, without any plasma particles.\n");
    print0("------------------------------------------------------------\n");
#endif

#if USE_PARTICLE_TRACING
#if !NTRL_ONLY
    //if Neutral only run, default is to trace neutrals which is usually expensive
    //default is electrons in global_grid
    particle_tracing::init_particle_tracing();
#endif
#endif


    /*************************************************************************
      Main Loop
     *************************************************************************/
    print0("\n            Beginning main loop  \n");
    auto begin_main_time = Benchmark::start();

    for (nstep = minstep; nstep < maxstep; ++nstep) 
    {
	printst0("---> nstep = %i \n",nstep);
		
#if USE_FEM_SOLVER
	printf(">> area weighted density calculation at %i ...\n", mpi_rank);
		area_density();
	printf(">> current through cell faces at %i ...\n", mpi_rank);
		cellfacecurrent();
	printf(" done!\n");
#endif

        if (is_statistics_step()){
	  global_grid.print_size_all();
        }

	particle_tracing::process(nstep);

	if(is_domain_recreation_step()){
	  print0("re-creating domain decomposition at nstep=%d \n",nstep);
	  init_domain_decomposition();
	  decompose_domain_all(global_grid);
	}
#if !NTRL_ONLY
#if GEOM_TYPE_HEMP
	MAX_dens_accumulate(global_grid.layers[ELECTRONS], MAX_Ne_NRmax, MAX_Ne_NZmax, &MAX_Ne_acc, &MAX_Ne_cntr);
#endif
#endif
#if DEBUG_PRTL_POS
	// TODO does not work we dont have any species anymore 
	//for (int i = 0; i < MAX_SPECIES; ++i){
	//  print_domain( nstep, array_labels[i], species[i].particles, *species[i].number_of_particles, mpi_rank );
	//}
	static_assert( false && "not implemented" );
#endif
	/*************************************************************************
	  Pushers
	 *************************************************************************/
	// always push electrons
	push_electrons(&pushE_time, Ampl);

	// push ions in ion steps
	if (is_ion_step()) {
	  push_ions( );
	}

	// push neutrals in neutral steps
#if !NTRL_CONST
	if (is_ntrl_step()){ 
	  push_neutrals( &pushN_time, &injN_time );
	}
#endif
	/*************************************************************************
	  Neutral Collisions
	 *************************************************************************/
#if !NTRL_CONST
	if (ncoll_n_n > 0 && nstep / ncoll_n_n * ncoll_n_n == nstep){
	  printst0("ntrl collisions\n");
	  //just a hotfix to make it work, will be deleted later
	  global_grid.layers[NEUTRALS].order();
	  neutral_collisions( 
			      &engcheck, 
			      &momcheck, 
			      &colN_time,
			      &orderN_time,
			      neutral_reaction_type
	  );
	  printst0("ntrl collisions finished\n");
	}
#endif
	/*************************************************************************
	  Collisions
	 *************************************************************************/
#if !NTRL_ONLY
	if (ncoll_n > 0 && nstep / ncoll_n * ncoll_n == nstep) {
	  printst0("charged particle collisions\n");
	  //just a hotfix to make it work, will be deleted later
	  global_grid.order_all();
	  collisions(  
		      &engcheck, 
		      &momcheck, 
		      &coll_time,
		      &order_time
	  );
	  printst0("charged particle collisions finished\n");
	}
#endif
	/*************************************************************************
	  Potential calculation
	 *************************************************************************/
#if !NTRL_ONLY
	printst0("potential calculation\n");
	calculate_potential( &field_time, r1,r2,z1,z2 ); 
	electric_field(phi, E_grid, E_ion);
	printst0("field calculation finished\n");
#endif
        /*************************************************************
		Injection
	*************************************************************/
	injection(&inj_time);

	/*     calculation and output the averaged values     */
	if (nstep > nav_start) {
#if !NTRL_ONLY
	  /* Phase Resolved */
	  //TODO since this is only relevant for RF move it somewhere else
#if GEOM_TYPE_RF
#if DIAG_PR
	    print0("collect data for phase resolved diagnostics\n");
	    global_grid.for_all_layers( [](GridLayer& layer){ 
	      if(layer.specie == ELECTRONS) {
		pr_coll(nstep, layer);
	      }
	      pr_diagnostics(nstep, layer);
	    } );
#endif
#endif
	  /**************************************************************************/
	  aver_moments(global_grid.layers[ELECTRONS]);
	  average(phi, phi_av );
#endif
	  if (nstep == nav_start+1) {
	    print0("Start averaging .......................nstep=%i \n",nstep);
	    print0(">>>Omega_pe*t= %-6.0f (%d steps) \n", nstep * dt, nstep);
#if !NTRL_ONLY
	    energy_all(global_grid.layers[ELECTRONS], v_te, En_e, En_f, E_grid, global_grid.mesh_dim, Ncell1, dt, dr);
	    print0("energy claculation done nstep=%i nav_start=%i nav_time=%i\n",nstep,nav_start,nav_time);
#endif
	  }
  
	  /****************************************************************
		            Output Diagnostics
	   ***************************************************************/
	  if (nstep == nav_start + nav_time)  
	  { 
	    print0("starting output ...\n");
	    write_diagnostics_cell_based();
	    reduce_for_output();
	    // CRITICAL DONT REDUCE phi_av
	    //reduce_double( phi_av, NG );
	    if ( mpi_rank == 0) {
	      write_diagnostics( 
		    field_time, 
		    pushE_time, 
		    coll_time, 
		    injN_time, 
		    pushN_time,
		    colN_time,
		    inj_time,
		    engcheck, 
		    begin_main_time,
		    momcheck,
		    program_options.datatrc
	     );
	   }
	   //print0("ni %i, nin %i\n", ni, nin );
	   //TODO here halt at diag via handler
           allreduce(&diag_halted);
           if (diag_halted) {
	     terminate_execution(  field_time, 
	       pushE_time, 
	       coll_time, 
	       injN_time,
	       pushN_time,
	       colN_time,
	       inj_time,
	       order_time,
	       orderN_time,
	       begin_main_time,
	       engcheck,
	       momcheck
	     );
	   }

           // all processes have to reset 
	   nav_start += nav_dt;
	   global_grid.for_all_layers( [](GridLayer& layer){ 
	     layer.n_aver = 0;
	   });
	   decompose_domain_all( global_grid );
	   /* -----  not DIAG  ----- */
	  }
	}

	/*********************************************
		sys_stop und back_up check
	**********************************************/
       if (nstep % BACKUP_EACH == 0 ) {
#if GEOM_TYPE_HEMP
#if !NTRL_ONLY
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      //TEMPORARY HOTFIX FOR SURFACE DENSITY SAVE
      //DELETE ONCE CELL SAVE STRUCTURE IS INTRODUCED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      for(int i=0; i<global_grid.z_dim;i++){
	surf_densZ[i]=0;
      }
      for(int specie=0; specie<MAX_SPECIES;specie++){
	GridLayer& layer=global_grid.layers[specie];
	for(int j=0; j<global_grid.z_dim;j++){
	  Cell& cell=layer.get_cell(nROUTthr, j);
	  surf_densZ[j]+=cell.surf_dens[1];
	  surf_densZ[j+1]+=cell.surf_dens[2];
	}
      }
      reduce(surf_densZ, global_grid.z_dim+1);
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif 
#endif 

#if GEOM_TYPE_EADS
#if !NTRL_ONLY
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      //TEMPORARY HOTFIX FOR SURFACE DENSITY READ FROM OLD RUNS TO MAP ONTO CELLS
      //DELETE ONCE CELL SAVE STRUCTURE IS INTRODUCED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      for(int i=0; i<global_grid.z_dim;i++){
	surf_densZ[i]=0;
      }
      for(int i=0; i<global_grid.r_dim;i++){
	surf_densR[i]=0;
      }
      std::cout<<"Adding loaded surface density charge in cells at r="<<nROUTthr<<std::endl;
      for(int specie=0; specie<MAX_SPECIES;specie++){
	GridLayer& layer=global_grid.layers[specie];
	for(int j=0; j<global_grid.z_dim;j++){
	  Cell& cell=layer.get_cell(nrMAXmodel, j);
	  surf_densZ[j]+=cell.surf_dens[1];
	  surf_densZ[j+1]+=cell.surf_dens[2];
	}
	for(int j=0; j<global_grid.r_dim;j++){
	  Cell& cell=layer.get_cell(j, nzMINmodel);
	  surf_densR[j]+=cell.surf_dens[1];
	  surf_densR[j+1]+=cell.surf_dens[0];
	}
      }
      reduce(surf_densZ, global_grid.z_dim+1);
      reduce(surf_densR, global_grid.r_dim+1);
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif 
#endif 
          reduce_for_output();
	  switch_backup_file = !switch_backup_file;
	  std::string filename = "save/particle_backup" + std::to_string(switch_backup_file) + ".h5";
	  h5_backup_particle_data(filename.c_str());
	  decompose_domain_all( global_grid );
       }


    //Add subcycling step
    if(nstep==ntrlstep){
      ntrlstep += global_grid.layers[NEUTRALS].dt_subcycling;
    }
#if !NTRL_ONLY
    if(nstep==ionstep){
      ionstep += global_grid.layers[POSITIVE_IONS].dt_subcycling;
    }
#endif

    } /* End of main loop  */

    nstep=maxstep-1 ;
    reduce_for_output();
    terminate_execution(  field_time, 
			  pushE_time, 
			  coll_time, 
			  injN_time,
			  pushN_time,
			  colN_time,
			  inj_time,
			  order_time,
			  orderN_time,
			  begin_main_time,
			  engcheck,
			  momcheck
    );
}
