/***************************************************************************

 Here we collect all output and diagnostic routines that are cell based
  
***************************************************************************/
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <string>
#include <iostream>
#include <fstream>

#include "var_dim.h"
#include "density.h"
#include "var_diag.h"
#include "var.h"
#include "pic.h"
#include "arrays.h"
#include "pic.h"
#include "grid.h"
#include "init_diag.h"
#include "macros.h"
#include "output.h"
#include "fem_solver.h"
#include "engy.h"
#include "emission.h"
#include "mpi_wrapper.h"
#include "rng64.h"
#include "emission.h"
#include "angular_current_diag.h"

#define DEBUG_LEVEL DEBUG_ERROR //DEBUG_ERROR  DEBUG_INFO_1
#include "debug_printing.h"

using namespace std;

void write_diagnostics(
    double field_time, 
    double pushE_time, 
    double coll_time, 
    double injN_time, 
    double pushN_time,
    double colN_time,
    double inj_time,
    double engcheck, 
    Bench_Start begin_main_time,
    Vec3d momcheck,
    std::string datatrc
)
{
  const int MAX_CHAR_SIZE = 256;
  char numb[MAX_CHAR_SIZE];

  snprintf(numb, MAX_CHAR_SIZE, "%08i", nstep);

  // the folder to which all output files will be stored
  std::string path = "out/";

  // this lambda will create string that represents the full path to the file and append ".dat" to it
  auto get_full_path = [&numb,&path]( std::string filename ){ 
    // this variable needs to be static so that we can return a pointer
    // and dont care about the memory that we otherwise would have to delete 
    std::string storage;
    storage = path + filename + std::string(numb) + ".dat";
    return storage; 
  };

  std::cout<<"Start global output"<<std::endl;

#if !NTRL_ONLY
  out_and_zero_phi(phi_av, dt, dr, get_full_path("phi")); 
#endif


#if !NTRL_ONLY
  //TODO add similar diagnostics for oxygen @paulm
  out_and_zero_coll(ncoll_coulomb, get_full_path("coll_coulomb") );
  out_and_zero_coll(ncoll_ioniz, get_full_path("coll_ioniz"));
  out_and_zero_coll(ncoll_el_ntrl, get_full_path("coll_el_ntrl"));
  out_and_zero_coll(ncoll_ion_ntrl, get_full_path("coll_ion_ntrl"));
  out_and_zero_coll(ncoll_el_ntrl_exc, get_full_path("coll_el_ntrl_exc"));
#if USE_FEEDGAS_OXYGEN
  out_and_zero_coll(ncoll_nion_ntrl, get_full_path("coll_nion_ntrl"));
  out_and_zero_coll(ncoll_dat, get_full_path("coll_dat"));
  out_and_zero_coll(ncoll_recomb, get_full_path("coll_recomb"));
  out_and_zero_coll(ncoll_ntrlz, get_full_path("coll_ntrlz"));
  out_and_zero_coll(ncoll_ntrl_detach, get_full_path("coll_ntrl_detach"));
  out_and_zero_coll(ncoll_el_det, get_full_path("coll_el_det"));
#endif
#if USE_TWOPLUS_IONS
  out_and_zero_coll(ncoll_double_ioniz, get_full_path("coll_double_ioniz"));
  out_and_zero_coll(ncoll_ioniz2, get_full_path("coll_ion_ioniz"));
#endif
#endif
#if !NTRL_CONST
  out_and_zero_coll(ncoll_ntrl_ntrl, get_full_path( "coll_ntrl_ntrl") );
#endif
  // TODO: adapt this vor varaible number of species (ions2p and negativ ions are missing!!)
  out_mfp(ecoll_dr,ecoll_dt,icoll_dr,icoll_dt,ncoll_dr,ncoll_dt, get_full_path("coll_mfp") ,n_bins_coll);


#if DIAG_ANGULAR_CURRENT
  for(auto i=0;i<MAX_SPECIES;++i){
    int k = POSITIVE_IONS;
#if USE_TWOPLUS_IONS
    k = TWOPLUS_IONS;
#endif
    if(i==POSITIVE_IONS || i==k){
      GridLayer& layer = global_grid.layers[i];
      std::string specie = layer.out_name;
      double nav_steps = static_cast<double>(nav_time);
      layer.angular_current.print_all(get_full_path(specie+"_angular_current_"),
					get_full_path(specie+"_angular_energy_"),
					get_full_path(specie+"_angular_location_"),nav_steps);
    }
  }
#endif

#if !NTRL_ONLY
#if !USE_HDF5
  energy_all(electrons, ne, v_te, &En_e, &En_f, E_grid, NG, Ncell1, dt, dr);
#endif
#endif

  // Trace output
  FILE *file_out;
  file_out = fopen(datatrc.c_str(), "a");
#if GEOM_TYPE_HEMP
#if !NTRL_ONLY
  out_and_zero_current_to_walls(count_particles);

#if USE_TWOPLUS_IONS
  fprintf(file_out, "%13.5e %8d %8d %8d %8d %13.5e %8d \n",
  		  dt_0 * nstep, nstep, ne_tmp, ni_tmp, nn_tmp, MAX_Ne, ni2p_tmp);
#else
  fprintf(file_out, "%13.5e %8d %8d %8d %8d %13.5e \n",
  		  dt_0 * nstep, nstep, ne_tmp, ni_tmp, nn_tmp, MAX_Ne);
#endif
#else 
  fprintf(file_out, "%13.5e %8d %8d %8d %8d %13.5e \n",
  		  dt_0 * nstep, nstep, ne_tmp, ni_tmp, nn_tmp, MAX_Ne );
#endif
#else
#if USE_FEEDGAS_OXYGEN
  fprintf(file_out, "%13.5e %8d %8d %8d %8d %8d %13.5e \n", dt_0 * nstep, nstep,
	  ne_tmp, ni_tmp, nin_tmp, nn_tmp, MAX_Ne);
#elif USE_TWOPLUS_IONS
  fprintf(file_out, "%13.5e %8d %8d %8d %8d %13.5e %8d \n", dt_0 * nstep,
	  nstep, ne_tmp, ni_tmp, nn_tmp, MAX_Ne, ni2p_tmp );
#else
  fprintf(file_out, "%13.5e %8d %8d %8d %8d %13.5e \n", dt_0 * nstep, nstep,
	  ne_tmp, ni_tmp, nn_tmp, MAX_Ne);
#endif
#endif
  fclose(file_out);

  if( get_rank() == 0) {
    //use ne_tmp, etc variables in trace output as this ensures the right particle number in any case
    printf(">>>Omega_pe*t= %-6.0f (%d steps) \n", nstep*dt, nstep);
    printf("ne = %d  ni = %d  nin = %d  ni2p = %d  nn = %d \n",ne_tmp,ni_tmp,nin_tmp,ni2p_tmp,nn_tmp);
    printf("............................................. End averaging \n");
    
    printf("Field: %11.3e       PushE: %11.3e       Collisions: "
     "%11.3e   InjN: %11.3e       PushN: %11.3e       ColN: "
     "%11.3e   Inject: %11.3e [ms]\n",field_time, pushE_time, coll_time, injN_time, pushN_time,
     colN_time, inj_time);
  	
	ofstream output( "benchmark.csv",  ios::app );
	output << nstep << "," << field_time << "," <<  pushE_time << "," <<  coll_time << ","; 
	output << injN_time << "," << pushN_time << "," <<  colN_time << "," << inj_time << endl;

    printf("\n");
    printf("EngCheck = %e  MomCheck.r = %e  MomCheck.t = %e  MomCheck.z = %e\n \n",
     engcheck, momcheck.r, momcheck.t, momcheck.z);

    double main_time = Benchmark::stop(begin_main_time);
    printf("Main_Time  %e ms \n", main_time );   
  }

  fflush(stdout);
}


void write_diagnostics_cell_based()
{
  const int MAX_CHAR_SIZE = 256;
  char numb[MAX_CHAR_SIZE];

  snprintf(numb, MAX_CHAR_SIZE, "%08i", nstep);

  // the folder to which all output files will be stored
  std::string path = "out/";

  // this lambda will create string that represents the full path to the file and append ".dat" to it
  auto get_full_path = [&numb,&path]( std::string filename ){ 
    // this variable needs to be static so that we can return a pointer
    // and dont care about the memory that we otherwise would have to delete 
    std::string storage;
    storage = path + filename + std::string(numb) + ".dat";
    return storage; 
  };

  for (auto& layer: global_grid.layers){
    string specie = layer.out_name;

    out_surf_dens(layer, get_full_path(specie+"_surf_dens_"));
    out_pt_too_fast(layer, get_full_path(specie+"_fast_pt_"));
//    out_charge_dens(layer, get_full_path(specie+"_charge_dens_"));

//density output
    out_dens(layer, get_full_path(specie+"_dens_"));

    //particle average velocities per cell  
    out_velz(layer, get_full_path( specie+"_velr_" ), get_full_path( specie+"_velt_" ), get_full_path( specie+"_velz_" ) );
    
    //output of Diagnostic velocity distribution, makes r and z energy cuts 
    out_vel_fluidlines(layer, number_fluidlines, get_full_path(specie+"_fluidlines_")); 
    
    //particle flux diagnostic 
//    out_flux(layer, get_full_path(specie+"_fluidz_prtcls_"), get_full_path(specie+"_fluidr_prtcls_"), get_full_path(specie+"_fluidz_enrgy_"), get_full_path(specie+"_fluidr_enrgy"));

    //cout<<"output of mean energy and temperature"<<endl;
    out_tempz(layer, get_full_path(specie+"_temp_"), get_full_path(specie+"_temp_therm_"));

    out_dist(layer, get_full_path(specie+"_distzall_"), get_full_path(specie+"_distrall_"), get_full_path(specie+"_distcell_"));

    //Phase resolved output routines (just important for RF)
    out_dens_pr(layer, get_full_path(specie+"_pr_dens_"));
    out_velz_pr(layer, get_full_path(specie+"_pr_velr_"), get_full_path(specie+"_pr_velt_"), get_full_path(specie+"_pr_velz_")); // TODO bug
    out_and_zero_dist_pr(layer, get_full_path(specie+"_pr_distzall_"), get_full_path(specie+"_pr_distrall_"), get_full_path(specie+"_pr_distcell_"));
  }

#if USE_FEM_SOLVER
		out_area_density(get_full_path("aw_charge_"));
		out_cellcurrent(get_full_path("cell_current_"));
    out_efield(get_full_path("fem_efield_"), get_full_path("stand_efield_"));
#endif

#if DIAG_EMISSION
  //Emission diagnostic, mainly used in eads
   emission("./sigma.dat", global_grid.layers[ELECTRONS], global_grid.layers[NEUTRALS], 
    get_full_path("emission"), 1e40);
#endif

#if !NTRL_ONLY
#if DIAG_COLL
  //this is just a test_case scenario how a layer based collision diag could look like
  out_coll_pr(global_grid.layers[ELECTRONS], get_full_path("ncoll_ionizpr"));
#endif
#endif

#if DIAG_PR
//reset all pr diag arrays, check if this is the easiest way to do it
//TODO add recent added cell diagnostic properties
for (int kind=0; kind<MAX_SPECIES; kind++){
  for(int i =0; i<global_grid.layers[kind].r_dim; i++){
    for(int j =0; j<global_grid.layers[kind].z_dim; j++){
      Cell& cell=global_grid.layers[kind].get_cell(i , j);
      for(int phase=0; phase<cell.diagnostic_arrays.pr_count.size();phase++){
	cell.diagnostic_arrays.particle_number_pr[phase]=0;
	cell.diagnostic_arrays.pr_count[phase]=0;
#if DIAG_VELZ
	cell.diagnostic_arrays.velzr_pr[phase]=0.;
	cell.diagnostic_arrays.velzz_pr[phase]=0.;
	cell.diagnostic_arrays.velzt_pr[phase]=0.;
#endif
#if DIAG_DIST
	for(int bin=0;bin<NBIN;bin++){
	  cell.diagnostic_arrays.vel_r_dist_pr[bin][phase]=0;
	  cell.diagnostic_arrays.vel_z_dist_pr[bin][phase]=0;
	  cell.diagnostic_arrays.vel_t_dist_pr[bin][phase]=0;
	  cell.diagnostic_arrays.vel_dist_pr[bin][phase]=0;
	}
#endif
      }
    }
  }
}
#endif


  //output of flux on a target (note: neglects negative ions and double positive ions)
  //+ output of energy flux to the model
#if !NTRL_ONLY
#if !NTRL_CONST
//  out_flux_P_sum(global_grid.layers[ELECTRONS], global_grid.layers[POSITIVE_IONS], global_grid.layers[NEUTRALS], nrMAXmodel,nzMINmodel, get_full_path( "Pflux_sum" ));
//  out_flux_P_sum(global_grid.layers[ELECTRONS], global_grid.layers[POSITIVE_IONS], global_grid.layers[NEUTRALS], nrMAXmodel,pos_flux_diag, get_full_path( "Pflux1_sum" ));
//  out_flux_E_sum(global_grid.layers[ELECTRONS], global_grid.layers[POSITIVE_IONS], global_grid.layers[NEUTRALS], nrMAXmodel,nzMINmodel, get_full_path( "Eflux_sum" ));
//  out_flux_E_sum(global_grid.layers[ELECTRONS], global_grid.layers[POSITIVE_IONS], global_grid.layers[NEUTRALS], nrMAXmodel,pos_flux_diag, get_full_path( "Eflux1_sum" ));
//#endif
//  out_flux_P(global_grid.layers[ELECTRONS], nrMAXmodel, nzMINmodel, get_full_path( "Pflux_e" ));
//  out_flux_P(global_grid.layers[ELECTRONS], nrMAXmodel, pos_flux_diag, get_full_path( "Pflux1_e" ));
//  
//  out_flux_P(global_grid.layers[POSITIVE_IONS], nrMAXmodel, nzMINmodel, get_full_path( "Pflux_i" ));
//  out_flux_P(global_grid.layers[POSITIVE_IONS], nrMAXmodel, pos_flux_diag, get_full_path( "Pflux1_i" ));
//  
//  out_flux_E(global_grid.layers[ELECTRONS], nrMAXmodel, nzMINmodel, get_full_path( "Eflux_e" ));
//  out_flux_E(global_grid.layers[ELECTRONS], nrMAXmodel, pos_flux_diag, get_full_path( "Eflux1_e" ));
//  
//  out_flux_E(global_grid.layers[POSITIVE_IONS], nrMAXmodel, nzMINmodel, get_full_path( "Eflux_i" ));
//  out_flux_E(global_grid.layers[POSITIVE_IONS], nrMAXmodel, pos_flux_diag, get_full_path( "Eflux1_i" ));
#endif
//
//#if !NTRL_CONST
//  out_flux_P(global_grid.layers[NEUTRALS], nrMAXmodel, nzMINmodel, get_full_path( "Pflux_n" ));
//  out_flux_P(global_grid.layers[NEUTRALS], nrMAXmodel, pos_flux_diag, get_full_path( "Pflux1_n" ));
//  out_flux_E(global_grid.layers[NEUTRALS], nrMAXmodel, nzMINmodel, get_full_path( "Eflux_n" ));
//  out_flux_E(global_grid.layers[NEUTRALS], nrMAXmodel, pos_flux_diag, get_full_path( "Eflux1_n" ));
#endif

  std::cout<<"Finished cell based output"<<endl;
}

void reduce_for_output(){
  printf("reducing all data to process 0\n");
#if !USE_HDF5
#if !NTRL_ONLY
  compose_domain( electrons, &ne, NE, NR, NZ); 
  compose_domain( ions, &ni, NE, NR, NZ);      //      it crashes for large number of particles (HEMP runs)
#endif
  compose_domain( neutrals, &nn, NN, NR, NZ);
#endif

#if USE_FEEDGAS_OXYGEN
#if !USE_HDF5
  compose_domain( nions, &nin, NE, NR, NZ);
#endif
#endif
#if USE_TWOPLUS_IONS
#if !USE_HDF5
  compose_domain( ions2p, &ni2p, NE, NR, NZ);
#endif
#endif

  //use temporary variables to store overall particle number
  //since in HDF5 particle numbers are not reduced 
#if !NTRL_ONLY
  ne_tmp=global_grid.layers[ELECTRONS].size();
  ni_tmp=global_grid.layers[POSITIVE_IONS].size();
#if USE_FEEDGAS_OXYGEN
  nin_tmp=global_grid.layers[NEGATIVE_IONS].size();
#endif
#if USE_TWOPLUS_IONS
  ni2p_tmp=global_grid.layers[TWOPLUS_IONS].size();
#endif
#endif
  nn_tmp=global_grid.layers[NEUTRALS].size();

#if USE_HDF5
#if USE_MPI
#if !NTRL_ONLY
  reduce(&ne_tmp);
  reduce(&ni_tmp);
#if USE_FEEDGAS_OXYGEN
  reduce(&nin_tmp);
#endif
#if USE_TWOPLUS_IONS
  reduce(&ni2p_tmp);
#endif
#endif
  reduce(&nn_tmp);


  for( auto& element : count_particles){
    reduce(&element.anode);  
    reduce(&element.dielectric);  
    reduce(&element.dielectric_SEE);  
    reduce(&element.metal);  
    reduce(&element.domain_boundary);  
  }

#if DIAG_ANGULAR_CURRENT
  for(auto i=0;i<MAX_SPECIES;++i){ 
    GridLayer& layer = global_grid.layers[i];
    layer.angular_current.reduce_for_output();
  }
#endif

#endif//USE_MPI
#endif//USE_HDF5
}


void out_dens(GridLayer& layer, std::string dat_name){
  if (layer.n_aver<1) layer.n_aver = 1;

  std::vector<double>  number_of_particles(layer.r_dim*layer.z_dim);

  for ( unsigned int r=0; r<layer.r_dim;r++){
      for ( unsigned int z=0; z<layer.z_dim;z++){
	  int ig=r*layer.z_dim+z;
	  Cell& cell=layer.get_cell(r , z);
	  number_of_particles[ig]=cell.diagnostic_arrays.particle_number; 
      }
  }
  
  reduce(number_of_particles);
  
  int mpi_rank=get_rank();
  if (mpi_rank==0){
    FILE  *file = fopen(dat_name.c_str(), "w");
    if ( !file ){ 
      printf("could not open %s\n",dat_name.c_str()); 
      exit(EXIT_FAILURE);
    }

    for ( unsigned int r = 0; r < layer.r_dim; ++r){
       double fi = 1./(Ncell1*(r+0.5)*layer.n_aver);    // = old version (due to NGP) 
       for (unsigned int z = 0; z < layer.z_dim; ++z){
          int ig=r*layer.z_dim+z;
          fprintf(file, "% .5e ",  fi*number_of_particles[ig]);               
          }
       fprintf(file, "\n");        
    }
    fclose(file);
  }
}

void out_dens_pr(GridLayer& layer, string dat1){
#if DIAG_PR

  Cell& init_cell = layer.get_cell(0,0); //obtain number of phases
  for (int phase=0; phase<init_cell.diagnostic_arrays.particle_number_pr.size(); phase++){
    //cout<<"Pr Diagnostic for phase "<<phase<<endl;
    std::vector<int>  number_of_particles(layer.r_dim*layer.z_dim);
    std::vector<int>  number_of_counts(layer.r_dim*layer.z_dim);

    for (int r=0; r<layer.r_dim;r++){
      for (int z=0; z<layer.z_dim;z++){
	int ig=r*layer.z_dim+z;
	Cell& cell=layer.get_cell(r , z);
	number_of_particles[ig]=cell.diagnostic_arrays.particle_number_pr[phase]; 
	number_of_counts[ig]=cell.diagnostic_arrays.pr_count[phase]; 
      }
    }
    reduce(number_of_particles);
    reduce(number_of_counts);
    
    int mpi_rank=get_rank();
    if (mpi_rank==0){
      FILE *f;
      std::string dat="out/" + dat1 + to_string(phase) + ".dat";
      f = fopen(dat.c_str(), "w");

      for (int r = 0; r < layer.r_dim; ++r){
	 for (int z = 0; z < layer.z_dim; ++z){
	    Cell& cell = layer.get_cell(r, z);
	    int ig=r*layer.z_dim+z;
	    double fi;
	    if (number_of_counts[ig]==0) {
	      fi=0;
	    }else{
	      fi = 1./(Ncell1*(r+0.5)*number_of_counts[ig]); 
	    }
	    fprintf(f, "% .5e ",  fi*number_of_particles[ig]);               
	    }
	 fprintf(f, "\n");        
      }
     fclose(f);
    }

  }
#endif
}


void out_velz(GridLayer& layer, string dat1, string dat2, string dat3) {
#if DIAG_VELZ
    double u0=layer.cs;
    u0 = MAX(u0, 1.e-10);

    std::vector<int> pn(layer.r_dim * layer.z_dim);
    std::vector<double> vel_r(layer.r_dim * layer.z_dim);
    std::vector<double> vel_t(layer.r_dim * layer.z_dim);
    std::vector<double> vel_z(layer.r_dim * layer.z_dim);

    for (int i=0; i<layer.r_dim; ++i){
      for (int j=0; j<layer.z_dim; ++j){
	Cell& cell = layer.get_cell(i , j);
	  int ind = i * layer.z_dim + j;
	  pn[ind]=cell.diagnostic_arrays.particle_number;
	  if (cell.diagnostic_arrays.particle_number>1){
	    vel_r[ind]=cell.diagnostic_arrays.particle_velocity[0];
	    vel_t[ind]=cell.diagnostic_arrays.particle_velocity[1];
	    vel_z[ind]=cell.diagnostic_arrays.particle_velocity[2];
	  }else{
	    vel_r[ind]=0.;
	    vel_t[ind]=0.;
	    vel_z[ind]=0.;
	  }
      }
    }
    reduce(pn);
    reduce(vel_r);
    reduce(vel_t);
    reduce(vel_z);

    int mpi_rank = get_rank();
    if(mpi_rank==0){
      FILE *file1, *file2, *file3;

	file1 = fopen(dat1.c_str(), "w");
	file2 = fopen(dat2.c_str(), "w");
	file3 = fopen(dat3.c_str(), "w");

	for (int i=0; i<layer.r_dim; i++){
	    for (int j=0; j<layer.z_dim; j++){
	      int ind = i * layer.z_dim + j;
	      if (pn[ind]>1){
		  vel_r[ind]/=(u0 * pn[ind]);
		  vel_t[ind]/=(u0 * pn[ind]);  
		  vel_z[ind]/=(u0 * pn[ind]);
		} else {
		  vel_r[ind]=vel_t[ind]=vel_z[ind]=0;
	      }
	      fprintf(file1, "% .5e ", vel_r[ind]); //in units of cs 
	      fprintf(file2, "% .5e ", vel_t[ind]); //in units of cs 
	      fprintf(file3, "% .5e ", vel_z[ind]); //in units of cs 
	    }
	    fprintf(file1, "\n");
	    fprintf(file2, "\n");
	    fprintf(file3, "\n");
	}

	fclose(file1);
	fclose(file2);
	fclose(file3);
    }
#endif
}

void out_velz_pr(GridLayer& layer, string datr, string datz, string datt) {

#if DIAG_PR
#if DIAG_VELZ
    double u0=layer.cs;
    u0 = MAX(u0, 1.e-10);

    Cell& init_cell = layer.get_cell(20,20); //TODO this is not very elegant, try to get the diag sizes from somewhere else
    for (int phase=0; phase<init_cell.diagnostic_arrays.velzr_pr.size(); phase++){
    std::vector<int> pn(layer.r_dim * layer.z_dim);
    std::vector<int> number_of_counts(layer.r_dim * layer.z_dim);
    std::vector<double> vel_r(layer.r_dim * layer.z_dim);
    std::vector<double> vel_t(layer.r_dim * layer.z_dim);
    std::vector<double> vel_z(layer.r_dim * layer.z_dim);

    for (int i=0; i<layer.r_dim; ++i){
      for (int j=0; j<layer.z_dim; ++j){
	Cell& cell = layer.get_cell(i , j);
	  int ind = i * layer.z_dim + j;
	  pn[ind]=cell.diagnostic_arrays.particle_number_pr[phase];
	  number_of_counts[ind]=cell.diagnostic_arrays.pr_count[phase]; 
	  if (cell.diagnostic_arrays.particle_number>1){
	    vel_r[ind]=cell.diagnostic_arrays.velzr_pr[phase];
	    vel_z[ind]=cell.diagnostic_arrays.velzz_pr[phase];
	    vel_t[ind]=cell.diagnostic_arrays.velzt_pr[phase];
	  }else{
	    vel_r[ind]=0.;
	    vel_z[ind]=0.;
	    vel_t[ind]=0.;
	  }
      }
    }
    reduce(pn);
    reduce(number_of_counts); 
    reduce(vel_r);
    reduce(vel_z);
    reduce(vel_t);

    int mpi_rank = get_rank();
    if(mpi_rank==0){

      std::string dat1="out/" + datr + to_string(phase) + ".dat";
      std::string dat2="out/" + datz + to_string(phase) + ".dat";
      std::string dat3="out/" + datt + to_string(phase) + ".dat";

	FILE *file1 = fopen(dat1.c_str(), "w");
	FILE *file2 = fopen(dat2.c_str(), "w");
	FILE *file3 = fopen(dat3.c_str(), "w");

	for (int i=0; i<layer.r_dim; i++){
	    for (int j=0; j<layer.z_dim; j++){
	      int ind = i * layer.z_dim + j;
	      if (pn[ind]>1){
		  vel_r[ind]/=(u0 * pn[ind]);
		  vel_z[ind]/=(u0 * pn[ind]);  
		  vel_t[ind]/=(u0 * pn[ind]);
		} else {
		  vel_r[ind]=vel_t[ind]=vel_z[ind]=0;
	      }
	      if(number_of_counts[ind]==0){ 
		fprintf(file1, "% .5e ", vel_r[ind]); //in units of cs 
		fprintf(file2, "% .5e ", vel_z[ind]); //in units of cs 
		fprintf(file3, "% .5e ", vel_t[ind]); //in units of cs 
	      } else {
		fprintf(file1, "% .5e ", vel_r[ind]/number_of_counts[ind]); //in units of cs 
		fprintf(file2, "% .5e ", vel_z[ind]/number_of_counts[ind]); //in units of cs 
		fprintf(file3, "% .5e ", vel_t[ind]/number_of_counts[ind]); //in units of cs 
	      }
	    }
	    fprintf(file1, "\n");
	    fprintf(file2, "\n");
	    fprintf(file3, "\n");
	}

	fclose(file1);
	fclose(file2);
	fclose(file3);
    }
  }
#endif
#endif
}



void out_tempz(GridLayer& layer, string dat1, string dat2) 
{
#if DIAG_TEMPZ
  double u0 = layer.cs;
  double fnorm=1.;
  if(layer.name=="electrons") fnorm = me_over_mi;

  FILE *file1, *file2;
  
  std::vector<int> pn(layer.r_dim*layer.z_dim);

  std::vector<double> vel_r(layer.r_dim*layer.z_dim);
  std::vector<double> vel_t(layer.r_dim*layer.z_dim);
  std::vector<double> vel_z(layer.r_dim*layer.z_dim);

  std::vector<double> temp_r(layer.r_dim*layer.z_dim);
  std::vector<double> temp_t(layer.r_dim*layer.z_dim);
  std::vector<double> temp_z(layer.r_dim*layer.z_dim);

  for (int i=0; i<layer.r_dim; ++i){
    for (int j=0; j<layer.z_dim; ++j){
      Cell& cell = layer.get_cell(i , j);
      int ind = layer.cell_number(i,j);
      auto& diag = cell.diagnostic_arrays;
      pn[ind]=diag.particle_number;

      temp_r[ind]=diag.particle_velocity_squared[0];
      temp_t[ind]=diag.particle_velocity_squared[1];
      temp_z[ind]=diag.particle_velocity_squared[2];

      vel_r[ind]=diag.particle_velocity[0];
      vel_t[ind]=diag.particle_velocity[1];
      vel_z[ind]=diag.particle_velocity[2];
    }
  }

  reduce(pn);
 
  reduce(vel_r);
  reduce(vel_t);
  reduce(vel_z);

  reduce(temp_r);
  reduce(temp_t);
  reduce(temp_z);

  int mpi_rank= get_rank();
  if (mpi_rank==0){
    for (int n = 0; n < layer.r_dim*layer.z_dim ; n++) {
      if (pn[n] >= 1) {
	temp_r[n] = temp_r[n]/(u0 * u0 * pn[n]) * fnorm;
	temp_t[n] = temp_t[n]/(u0 * u0 * pn[n]) * fnorm;
	temp_z[n] = temp_z[n]/(u0 * u0 * pn[n]) * fnorm;

	vel_r[n] = temp_r[n] - (SQU(vel_r[n]/(u0 * pn[n])) * fnorm);
	vel_t[n] = temp_t[n] - (SQU(vel_t[n]/(u0 * pn[n])) * fnorm);
	vel_z[n] = temp_z[n] - (SQU(vel_z[n]/(u0 * pn[n])) * fnorm);
      } else {
	vel_r[n] = vel_t[n] = vel_z[n] = 0.;
	temp_r[n] = temp_t[n] = temp_z[n] = 0.;
      }
    }
    file1 = fopen(dat1.c_str(), "w");
    file2 = fopen(dat2.c_str(), "w");
    for (int i = 0; i < layer.r_dim; ++i) {
      for (int j = 0; j < layer.z_dim; ++j) {
	fprintf(file1, "% .5e ", temp_r[i * layer.z_dim + j] + temp_t[i * layer.z_dim + j] + temp_z[i * layer.z_dim + j]);
	fprintf(file2, "% .5e ", vel_r[i * layer.z_dim + j] + vel_t[i * layer.z_dim + j] + vel_z[i * layer.z_dim + j]);
	//in units of Te_0
      }
      fprintf(file1, "\n");
      fprintf(file2, "\n");
    }
    fclose(file1);
    fclose(file2);
  }
#endif
}




void out_vel_fluidlines(GridLayer& layer, int n, string dat) {
#if DIAG_FLUID
#if DIAG_VELZ
  int tmax = 1000; //maximum number of possible steps in current domain

  FILE *fout;
  fout = fopen(dat.c_str(), "w");

  std::vector<int> pn(layer.r_dim*layer.z_dim);
  std::vector<double> velr(layer.r_dim*layer.z_dim);
  std::vector<double> velz(layer.r_dim*layer.z_dim);

  // unnormalized averaged velocity
  for (int i = 0; i < layer.r_dim; i++) {
    for (int j = 0; j < layer.z_dim; j++) {
      int ind = i*layer.z_dim + j;
      Cell& cell = layer.get_cell(i , j);
      pn[ind] = cell.diagnostic_arrays.particle_number;
      velr[ind] = cell.diagnostic_arrays.particle_velocity[0];
      velz[ind] = cell.diagnostic_arrays.particle_velocity[2];
    }
  }

  reduce(pn);
  reduce(velr);
  reduce(velz);

  int mpi_rank = get_rank();
  if(mpi_rank==0){
    
    double vr, vz;
    std::vector<double> r;
    std::vector<double> z;
    r.resize(n);
    z.resize(n);

    // integration of averaged velocities
    for (int t = 0; t < tmax; ++t) {
      for (int i = 0; i < n; ++i) {
	if (t == 0) {
	  r[i] = RAND * 50.; //region where trajectories start
	  z[i] = RAND * 50.;
	}

	int indi = (int)r[i];
	int indj = (int)z[i];

	if (indi >= layer.r_dim || indj >= layer.z_dim || indi < 0 || indj < 0) {
	  vr = vz = 0.;
	} else {
	  vr = velr[indi * layer.z_dim + indj];
	  vz = velz[indi * layer.z_dim + indj];

	  vr=MAX(vr,1e-12); 	
	  vz=MAX(vz,1e-12); 	
	  vr = vr / sqrt(vr * vr + vz * vz);
	  vz = vz / sqrt(vr * vr + vz * vz);
	}
	r[i] += vr;
	z[i] += vz;

	fprintf(fout, "%i %e %e\n", i, r[i], z[i]);
      }
    }

    fclose(fout);
  }
#endif
#endif
}




void out_flux(GridLayer& layer, string dat1, string dat2, string dat3, string dat4)
{ 
#if DIAG_FLUX
  int vector_size = (layer.r_dim+1) * (layer.z_dim+1); 
  std::vector<double>   fluidr(vector_size);
  std::vector<double>   fluidz(vector_size);
  std::vector<double>   energy_fluidr(vector_size);
  std::vector<double>   energy_fluidz(vector_size);

  for (int i = 0; i < layer.r_dim; ++i){
    for (int j = 0; j < layer.z_dim; ++j){   
      Cell& cell=layer.get_cell(i , j);
      fluidr[(i+1)*(layer.z_dim+1) + j] -= cell.diagnostic_arrays.particle_flux[0]; 
      fluidr[i*(layer.z_dim+1) + j] += cell.diagnostic_arrays.particle_flux[2]; 
      fluidz[i*(layer.z_dim+1) + j + 1] -= cell.diagnostic_arrays.particle_flux[1]; 
      fluidz[i*(layer.z_dim+1) + j] += cell.diagnostic_arrays.particle_flux[3]; 

      energy_fluidr[(i+1)*(layer.z_dim+1) + j] -= cell.diagnostic_arrays.energy_flux[0]; 
      energy_fluidr[i*(layer.z_dim+1) + j] += cell.diagnostic_arrays.energy_flux[2]; 
      energy_fluidz[i*(layer.z_dim+1)+ j + 1] -= cell.diagnostic_arrays.energy_flux[1]; 
      energy_fluidz[i*(layer.z_dim+1) + j] += cell.diagnostic_arrays.energy_flux[3]; 
    }
  }

  reduce(fluidr);
  reduce(fluidz);
  reduce(energy_fluidr);
  reduce(energy_fluidz);

  int mpi_rank=get_rank();
  if (mpi_rank==0){
    FILE    *fz, *fr, *er, *ez;
    fz = fopen(dat1.c_str(), "w");
    fr = fopen(dat2.c_str(), "w");
    ez = fopen(dat3.c_str(), "w");
    er = fopen(dat4.c_str(), "w");
    
    double m  = layer.mass*elemM;           // particle mass in [kg]
    double dt = dt_0*layer.dt_subcycling;   // subcycling time in [s]
    double dr = dr_0*1e-2;                  // cell size dr_0 in [m]
    iprintf(" m=%e (layer.mass=%f) dt=%f (layer.dt_subcycling=%i) dr=%f\n",m, layer.mass,dt,layer.dt_subcycling,dr);
    iprintf(" elemM=%e dt_0=%e dr_0=%e N_SP=%f nav_time=%i\n", elemM, dt_0,dr_0,N_SP,nav_time);

    for (int i = 0; i < layer.r_dim+1; ++i){
      //calculate area for scaling in r and z direction [m^2] 
      double r=i*dr;
      double R_in=r;
      double R_out=r+dr;
      double Az=PI*(SQU(R_out)-SQU(R_in)); 
      double Ar=TWOPI*(i*dr_0)*dr_0;   // = TWOPI * radius of zylinder * width of zylinder
       for (int j = 0; j < layer.z_dim+1; ++j){
          if (Ar!=0){
            fprintf(fr, "% .5e ",  fluidr[i*(layer.z_dim+1) + j]*N_SP/Ar/nav_time/dt_0);  // in [1/m^2/s]
            fprintf(er, "% .5e ",  energy_fluidr[i*(layer.z_dim+1) + j]*N_SP*SQU(dr_0/dt)*m*0.5/Ar/nav_time/dt_0); // in [J/m^2/s]
          }else{
            fprintf(fr, "% .5e ",  0.);               
            fprintf(er, "% .5e ",  0.);               
          }
          fprintf(fz, "% .5e ",  fluidz[i*(layer.z_dim+1) + j]*N_SP/Az/nav_time/dt_0);    // in [1/m^2/s] 
          fprintf(ez, "% .5e ",  energy_fluidz[i*(layer.z_dim+1) + j]*N_SP*SQU(dr_0/dt)*m*0.5/Az/nav_time/dt_0); // in [J/m^2*/s]
          }
       fprintf(fr, "\n");        
       fprintf(fz, "\n");        
       fprintf(er, "\n");        
       fprintf(ez, "\n");        
    }

    fclose(fr); 
    fclose(fz); 
    fclose(er); 
    fclose(ez); 
  }
#endif
}


void out_dist(GridLayer& layer, string fnDistz, string fnDistr, string fnDist)
{
#if DIAG_DIST
  FILE    *fr,*fz, *f;
  
  std::vector<int> distr_n(layer.z_dim);
  std::vector<int> distz_n(layer.r_dim);
  std::vector<int> dist_n_all(layer.r_dim*layer.z_dim);
  std::vector<vector<int>> dist_all(layer.r_dim*layer.z_dim);

  std::vector<vector<int>> distr_r(layer.z_dim);
  std::vector<vector<int>> distr_z(layer.z_dim);
  std::vector<vector<int>> distr_t(layer.z_dim);

  std::vector<vector<int>> distz_r(layer.r_dim);
  std::vector<vector<int>> distz_z(layer.r_dim);
  std::vector<vector<int>> distz_t(layer.r_dim);

  for (int j=0; j<layer.z_dim; ++j){
    distr_r[j].resize(NBIN);
    distr_z[j].resize(NBIN);
    distr_t[j].resize(NBIN);
  }

  for (int i=0; i<layer.r_dim; ++i){
    distz_r[i].resize(NBIN);
    distz_z[i].resize(NBIN);
    distz_t[i].resize(NBIN);
    for (int j=0; j<layer.z_dim; ++j){
      int ind = i * layer.z_dim + j;
      dist_all[ind].resize(NBIN);
    }
  }

  //initialize all arrays to zero
  for(int j=0; j<layer.z_dim; ++j){
    distr_n[j] = 0 ;
    for(int n=0; n<NBIN; ++n) {
      distr_r[j][n] = 0 ;
      distr_z[j][n] = 0 ;
      distr_t[j][n] = 0 ;
    }
  }
  for(int i=0; i<layer.r_dim; ++i){
    distz_n[i] = 0 ;
    for(int n=0; n<NBIN; ++n) {
      distz_r[i][n] = 0 ;
      distz_z[i][n] = 0 ;
      distz_t[i][n] = 0 ;
    }
    for(int j=0; j<layer.z_dim; ++j){
      int ind = i * layer.z_dim + j;
      dist_n_all[ind]= 0 ;
      for(int n=0; n<NBIN; ++n) {
	dist_all[ind][n] = 0 ;
      }
    }
  }
  //end initialize

  //cuts in r direction
  for(int j=0; j<layer.z_dim; ++j) {
    for(int i=0; i<layer.r_dim; ++i) {
      Cell& cell =layer.get_cell(i , j);
      distr_n[j] += cell.size();
      for(int n=0; n<NBIN;++n){
	if (cell.diagnostic_arrays.vel_r_dist.empty()){
	  distr_r[j][n]+=0;
	  distr_z[j][n]+=0;
	  distr_t[j][n]+=0;
	} else {
	  distr_r[j][n]+=cell.diagnostic_arrays.vel_r_dist[n];
	  distr_z[j][n]+=cell.diagnostic_arrays.vel_z_dist[n];
	  distr_t[j][n]+=cell.diagnostic_arrays.vel_t_dist[n];
	}
      }
    }
  }

  //cuts in r direction
  for(int i=0; i<layer.r_dim; ++i) {
    for(int j=0; j<layer.z_dim; ++j) {
      Cell& cell =layer.get_cell(i , j);
      distz_n[i] += cell.size();
      for(int n=0; n<NBIN;++n){
	if (cell.diagnostic_arrays.vel_r_dist.empty()){
	  distr_r[i][n]+=0;
	  distr_z[i][n]+=0;
	  distr_t[i][n]+=0;
	} else {
	  distz_r[i][n]+=cell.diagnostic_arrays.vel_r_dist[n];
	  distz_z[i][n]+=cell.diagnostic_arrays.vel_z_dist[n];
	  distz_t[i][n]+=cell.diagnostic_arrays.vel_t_dist[n];
	}
      }
    }
  }

  //overall distribution
  for(int i=0; i<layer.r_dim; ++i) {
    for(int j=0; j<layer.z_dim; ++j) {
      Cell& cell =layer.get_cell(i , j);
      int ind = i * layer.z_dim + j;
      dist_n_all[ind] = cell.size();
      for(int n=0; n<NBIN;++n){
	if (cell.diagnostic_arrays.vel_dist.empty()){
	  distr_r[i][n]+=0;
	} else {
	  dist_all[ind][n]=cell.diagnostic_arrays.vel_dist[n];
	}
      }
    }
  }

  reduce(distr_n);
  reduce(distz_n);

  reduce(distr_r);
  reduce(distr_z);
  reduce(distr_t);
  
  reduce(distz_r);
  reduce(distz_z);
  reduce(distz_t);

  reduce(dist_all);

  int mpi_rank = get_rank();
  if(mpi_rank==0){
    f = fopen(fnDist.c_str(), "w");
    fr = fopen(fnDistr.c_str(), "w");
    fz = fopen(fnDistz.c_str(), "w");
    if ( !fr || !fz || !f) {
      printf("could not open %s or %s or %s\n",fnDistr.c_str(),fnDistz.c_str(), fnDist.c_str());
      exit(EXIT_FAILURE);
    }

    for(int j=0; j<layer.z_dim; ++j) {
      for(int n=0;n<NBIN;++n){ 
	fprintf(fz, "%d %d %d %d %d %d \n",j,n,distr_n[j],
	    distr_r[j][n],distr_z[j][n],distr_t[j][n]);
      }
    }

    for(int i=0; i<layer.r_dim; ++i) {
      for(int n=0;n<NBIN;++n){ 
	fprintf(fr, "%d %d %d %d %d %d\n",i,n,distz_n[i],
	    distz_r[i][n],distz_z[i][n],distz_t[i][n]);
      }
    }

    for(int i=0; i<layer.r_dim; ++i) {
      for(int j=0; j<layer.z_dim; ++j) { 
	int ind = i * layer.z_dim + j;
	for(int n=0; n<NBIN; ++n) {
	  fprintf(f, "%d %d %d %d %d \n",i,j,n,dist_n_all[ind],dist_all[ind][n]);
	}
      }
    }

    fclose(f);
    fclose(fr);
    fclose(fz);
  }//mpi_rank==0

#endif
}

void out_and_zero_dist_pr(GridLayer& layer, string fnDistz, string fnDistr, string fnDist)
{
#if DIAG_PR
#if DIAG_DIST
  Cell& init_cell = layer.get_cell(20,20); //TODO this is not very elegant
  for (int phase=0; phase<init_cell.diagnostic_arrays.pr_count.size(); phase++){
  

  std::vector<int> distr_n(layer.z_dim,0);
  std::vector<int> distz_n(layer.r_dim,0);
  std::vector<int> dist_n_all(layer.r_dim*layer.z_dim,0);
  std::vector<vector<int>> dist_all(layer.r_dim*layer.z_dim);

  std::vector<vector<int>> distr_r(layer.z_dim);
  std::vector<vector<int>> distr_z(layer.z_dim);
  std::vector<vector<int>> distr_t(layer.z_dim);

  std::vector<vector<int>> distz_r(layer.r_dim);
  std::vector<vector<int>> distz_z(layer.r_dim);
  std::vector<vector<int>> distz_t(layer.r_dim);

  //initialize all arrays to zero
  for (int j=0; j<layer.z_dim; ++j){
    distr_r[j].resize(NBIN,0);
    distr_z[j].resize(NBIN,0);
    distr_t[j].resize(NBIN),0;
  }

  for (int i=0; i<layer.r_dim; ++i){
    distz_r[i].resize(NBIN,0);
    distz_z[i].resize(NBIN,0);
    distz_t[i].resize(NBIN,0);
    for (int j=0; j<layer.z_dim; ++j){
      dist_all[layer.cell_number(i,j)].resize(NBIN,0);
    }
  }
  //end initialize

  //cuts in r direction
  for(int j=0; j<layer.z_dim; ++j) {
    for(int i=0; i<layer.r_dim; ++i) {
      Cell& cell =layer.get_cell(i , j);
      distr_n[j] += cell.size();
      for(int n=0; n<NBIN;++n){
	distr_r[j][n]+=cell.diagnostic_arrays.vel_r_dist_pr[n][phase];
	distr_z[j][n]+=cell.diagnostic_arrays.vel_z_dist_pr[n][phase];
	distr_t[j][n]+=cell.diagnostic_arrays.vel_t_dist_pr[n][phase];
      }
    }
  }

  //cuts in z direction
  for(int i=0; i<layer.r_dim; ++i) {
    for(int j=0; j<layer.z_dim; ++j) {
      Cell& cell =layer.get_cell(i , j);
      distz_n[i] += cell.size();
      for(int n=0; n<NBIN;++n){
	distz_r[i][n]+=cell.diagnostic_arrays.vel_r_dist_pr[n][phase];
	distz_z[i][n]+=cell.diagnostic_arrays.vel_z_dist_pr[n][phase];
	distz_t[i][n]+=cell.diagnostic_arrays.vel_t_dist_pr[n][phase];
      }
    }
  }

  //overall distribution
  for(int i=0; i<layer.r_dim; ++i) {
    for(int j=0; j<layer.z_dim; ++j) {
      Cell& cell =layer.get_cell(i , j);
      int ind = i * layer.z_dim + j;
      dist_n_all[ind] = cell.size();
      for(int n=0; n<NBIN;++n){
	dist_all[ind][n]+=cell.diagnostic_arrays.vel_dist_pr[n][phase];
      }
    }
  }

  reduce(distr_n);
  reduce(distz_n);

  reduce(distr_r);
  reduce(distr_z);
  reduce(distr_t);

  reduce(distz_r);
  reduce(distz_z);
  reduce(distz_t);

  reduce(dist_all);

  int mpi_rank = get_rank();
  if(mpi_rank==0){
  FILE    *fr,*fz, *f;
    std::string dat1="out/" + fnDist  + to_string(phase) + ".dat";
    std::string dat2="out/" + fnDistr + to_string(phase) + ".dat";
    std::string dat3="out/" + fnDistz + to_string(phase) + ".dat";
    f = fopen(dat1.c_str(), "w");
    fr = fopen(dat2.c_str(), "w");
    fz = fopen(dat3.c_str(), "w");

    if ( !fr || !fz || !f) {
      printf("could not open %s or %s or %s\n",fnDistr.c_str(),fnDistz.c_str(), fnDist.c_str());
      exit(EXIT_FAILURE);
    }

    for(int j=0; j<layer.z_dim; ++j) {
      for(int n=0;n<NBIN;++n){ 
	fprintf(fz, "%d %d %d %d %d %d \n",j,n,distr_n[j],
	    distr_r[j][n],distr_z[j][n],distr_t[j][n]);
      }
    }

    for(int i=0; i<layer.r_dim; ++i) {
      for(int n=0;n<NBIN;++n){ 
	fprintf(fr, "%d %d %d %d %d %d\n",i,n,distz_n[i],
	    distz_r[i][n],distz_z[i][n],distz_t[i][n]);
      }
    }

    for(int i=0; i<layer.r_dim; ++i) {
      for(int j=0; j<layer.z_dim; ++j) { 
	int ind = i * layer.z_dim + j;
	for(int n=0; n<NBIN; ++n) {
	  fprintf(f, "%d %d %d %d %d \n",i,j,n,dist_n_all[ind],dist_all[ind][n]);
	}
      }
    }

    fclose(f);
    fclose(fr);
    fclose(fz);
  }//mpi_rank==0

  }//phases
#endif      
#endif      
}



void out_coll_pr(GridLayer& layer, string dat){
  
#if DIAG_PR
#if DIAG_COLL
  Cell& init_cell = layer.get_cell(20,20); //TODO this is not very elegant, try to get the diag sizes from somewhere else
  for (int phase=0; phase<init_cell.diagnostic_arrays.coll_ioniz_pr.size(); phase++){

    std::vector<double> coll(layer.r_dim*layer.z_dim);

    for (int i = 0; i < layer.r_dim; ++i) {
       for (int j = 0; j < layer.z_dim; ++j){
	  Cell& cell = layer.get_cell(i , j);
	  int ind =i*layer.z_dim +j;
	  if (cell.diagnostic_arrays.pr_count[phase]==0) {
	    coll[ind]=0.;
	  } else {
	    coll[ind]=(double) cell.diagnostic_arrays.coll_ioniz_pr[phase]/cell.diagnostic_arrays.pr_count[phase];
	  }
       }
   }

   reduce(coll);

   int mpi_rank=get_rank();
    if(mpi_rank==0){

      FILE    *f;
      std::string dat1="out/" + dat + to_string(phase) + ".dat";

      f = fopen(dat1.c_str(), "w");

      for (int i = 0; i < layer.r_dim; ++i){
	 for (int j = 0; j < layer.z_dim; ++j){
	    Cell& cell = layer.get_cell(i , j);
	    int ind=i*layer.z_dim +j;
	    fprintf(f, "%e ", coll[ind]);
	 }
	 fprintf(f, "\n");        
      }
      fclose(f); 
    }
  }
#endif
#endif
}

void out_flux_P(GridLayer& layer, int nr_model, int z_pos, string dat){
#if DIAG_FLUX  
  std::vector<double> flux(nr_model);  

   int mpi_rank=get_rank();
  for (int i=0; i<nr_model; i++){
    double r=i*dr_0;
    double R_in=r;
    double R_out=r+dr_0;
    double A=PI*(SQU(R_out)-SQU(R_in));
    int j = z_pos+1;
    Cell& cell = layer.get_cell(i,j);
    flux[i]=cell.diagnostic_arrays.particle_flux[3] / A * N_SP / dt_0 / nav_time ;
  }

  reduce(flux);

  if(mpi_rank==0){
    FILE     *f;
    f = fopen(dat.c_str(), "w");
    if ( f==NULL ) {printf("file not open %s \n",dat.c_str()); exit(-1);}
	 
    for (int i=0;i<nr_model;i++){
       fprintf(f, "%i  %20.16e \n",  i,flux[i]);  //particles per second for particles and W/s for energy              
    }
    fclose(f);
  }
#endif
}

void out_flux_E(GridLayer& layer, int nr_model, int z_pos, string dat){
#if DIAG_FLUX
  std::vector<double> flux(nr_model);  
  double m_e = elemM;
  double mass=layer.mass;
  double dt= dt_0*layer.dt_subcycling;
  int mpi_rank = get_rank();

  for (int i=0; i<nr_model; i++){
    double r=i*dr_0;
    double R_in=r;
    double R_out=r+dr_0;
    double A=PI*(SQU(R_out)-SQU(R_in));
    int j = z_pos+1;
    Cell& cell = layer.get_cell(i,j);
    flux[i]=cell.diagnostic_arrays.energy_flux[3] / A * N_SP / dt_0 / nav_time ;
    flux[i]*= SQU(dr_0/dt)*m_e*mass*0.5;
}

  reduce(flux);

  if(mpi_rank==0){
    FILE     *f;
    f = fopen(dat.c_str(), "w");
    if ( f==NULL ) {printf("file not open %s \n",dat.c_str()); exit(-1);}
	 
    for (int i=0;i<nr_model;i++){
       fprintf(f, "%i  %20.16e \n",  i,flux[i]);  //particles per second for particles and W/s for energy              
    }
    fclose(f);
  }
#endif
}

void out_flux_P_sum( GridLayer& el_layer, GridLayer& ion_layer, GridLayer& ntrl_layer, int nr_model, int z_pos, string dat_sum )
{
#if DIAG_FLUX
  std::vector<double> flux_el(nr_model);  
  std::vector<double> flux_ntrl(nr_model);  
  std::vector<double> flux_ion(nr_model);  
  std::vector<double> flux_sum(nr_model);  

    for (int i=0; i<nr_model; i++){
      double r=i*dr_0;
      double R_in=r;
      double R_out=r+dr_0;
      double A=PI*(SQU(R_out)-SQU(R_in));
      int j = z_pos+1;
      Cell& el_cell = el_layer.get_cell(i,j);
      Cell& ion_cell = ion_layer.get_cell(i,j);
      Cell& ntrl_cell = ntrl_layer.get_cell(i,j);
      flux_el[i]=el_cell.diagnostic_arrays.particle_flux[3]  ;
      flux_ion[i]=ion_cell.diagnostic_arrays.particle_flux[3] ;
      flux_ntrl[i]=ntrl_cell.diagnostic_arrays.particle_flux[3] * coll_fac_ntrl_ntrl ;
      flux_sum[i]= ( flux_el[i]+flux_ion[i]+flux_ntrl[i] ) / A * N_SP / dt_0 / nav_time;
    }

    reduce(flux_sum);

    int mpi_rank = get_rank();
    if(mpi_rank==0){
       FILE     *fsum;
       fsum = fopen(dat_sum.c_str(), "w");
       if ( fsum==NULL ) {printf("file not open %s \n",dat_sum.c_str()); exit(-1);}
	 
       for (int i=0;i<flux_sum.size();i++){
	 fprintf(fsum, "%i  %20.16e \n",  i,flux_sum[i] / nav_time);  //particles per second for particles and W/s for energy              
       }
       fclose(fsum);
    }
#endif
}

void out_flux_E_sum( GridLayer& el_layer, GridLayer& ion_layer, GridLayer& ntrl_layer, int nr_model, 
    int z_pos, string dat_sum )
{
#if DIAG_FLUX
  std::vector<double> flux_el(nr_model);  
  std::vector<double> flux_ntrl(nr_model);  
  std::vector<double> flux_ion(nr_model);  
  std::vector<double> flux_sum(nr_model);  

  double m_e = elemM; //mass electron in kg
  double dt_el=dt_0*el_layer.dt_subcycling;
  double dt_ion=dt_0*ion_layer.dt_subcycling;
  double dt_ntrl=dt_0*ntrl_layer.dt_subcycling;
 
    for (int i=0; i<nr_model; i++){
      double r=i*dr_0;
      double R_in=r;
      double R_out=r+dr_0;
      double A=PI*(SQU(R_out)-SQU(R_in));
      int j = z_pos+1;
      Cell& el_cell = el_layer.get_cell(i,j);
      Cell& ion_cell = ion_layer.get_cell(i,j);
      Cell& ntrl_cell = ntrl_layer.get_cell(i,j);
      flux_el[i]=el_cell.diagnostic_arrays.energy_flux[3] ;
      flux_ion[i]=ion_cell.diagnostic_arrays.energy_flux[3] ;
      flux_ntrl[i]=ntrl_cell.diagnostic_arrays.energy_flux[3] * coll_fac_ntrl_ntrl ;
      
      flux_el[i]*= SQU(dr_0/dt_el)*m_e;
      flux_ion[i]*= SQU(dr_0/dt_ion)*m_e*ion_layer.mass;
      flux_ntrl[i]*= SQU(dr_0/dt_ntrl)*m_e*ntrl_layer.mass;
      flux_sum[i]= 0.5 * ( flux_el[i]+flux_ion[i]+flux_ntrl[i]) / A * N_SP / dt_0 / nav_time;
    }

    reduce(flux_sum);

    int mpi_rank = get_rank();
    if(mpi_rank==0){
       FILE     *fsum;
       fsum = fopen(dat_sum.c_str(), "w");
       if ( fsum==NULL ) {printf("file not open %s \n",dat_sum.c_str()); exit(-1);}
	 
       for (int i=0;i<nr_model;i++){
	 fprintf(fsum, "%i  %20.16e \n",  i,flux_sum[i] / nav_time);  //particles per second for particles and W/s for energy              
       }
       fclose(fsum);
    }
#endif
}


//Non cell based diagnostic output functions
void out_and_zero_phi( double phi_av[], double dt, double dr, string dat1 )
{
  assert( nav_time>=1 && "nav_time = 0 steps");  
  double f1 = 2.*SQU(dr/dt)/nav_time*T_e0;     // remove 2. wenn you change qe and E 

  FILE    *file1 = fopen(dat1.c_str(), "w");
  for (int i = 0; i < global_grid.mesh_r_dim; ++i)  
  {
     for (int j = 0; j < global_grid.mesh_z_dim; ++j) 
        fprintf(file1, "% .5e ",  f1*phi_av[i*(global_grid.mesh_z_dim) + j]);  //volt             
     fprintf(file1, "\n");        
  }
  fclose(file1); 

  for (int i = 0; i < global_grid.mesh_r_dim; ++i)  {
    for (int j = 0; j < global_grid.mesh_z_dim; ++j) {
      phi_av[i*(global_grid.mesh_z_dim) + j] = 0.0;
    }
  }

}

//Non cell based diagnostic output functions
void out_rhs( double phi_av[], double dt, double dr, string dat1 )
{
  double f1 = 2.*SQU(dr/dt)/T_e0;     // remove 2. wenn you change qe and E 

  FILE    *file1 = fopen(dat1.c_str(), "w");
  for (int i = 0; i < global_grid.mesh_r_dim; ++i)  
  {
     for (int j = 0; j < global_grid.mesh_z_dim; ++j) 
        fprintf(file1, "% .5e ",  f1*phi_av[i*(global_grid.mesh_z_dim) + j]);  //volt             
     fprintf(file1, "\n");        
  }
  fclose(file1); 
}

void out_and_zero_coll( int ncoll[], string dat1 )
{
#if DIAG_COLL
  FILE    *file1;
  file1 = fopen(dat1.c_str(), "w");
  for (int i = 0; i <  global_grid.r_dim; ++i)  
  {
     for (int j = 0; j <  global_grid.z_dim; ++j) 
     {
        fprintf(file1, "%e ",  (double)ncoll[i*global_grid.z_dim + j] /nav_time);  //collisions per second             
        ncoll[i*global_grid.z_dim+j]=0;  //reset to zero after writing, to start new averaging timestep
     }
     fprintf(file1, "\n");        
  }
  fclose(file1); 
#endif
}


// Output routine for angular current diagnostic including energy distribution for each angle
void out_ang_curr_enDist( long int& ions_tot, std::vector<double>& Angle, std::vector<int>& Bin_particles, 
			  std::vector<double>& Energy,  std::vector<std::vector<int>>& Bin_particle_energy,
			  std::vector<std::vector<std::vector<int> >>& Bin_particle_origin,
			  double nav_steps, std::string f_name )
{
  double fact_pa_A = 1/(nav_steps*I_Amp_to_n2inj);  // scales from number of psudo particles -> Ampere
  double fact_pa = N_SP/(nav_steps);       // scales from number of psudo particles -> average number of particles / time step
  FILE    *file = fopen(f_name.c_str(), "w");
  if ( !file ){ 
    printf("could not open %s\n",f_name.c_str()); 
    exit(EXIT_FAILURE);
  }
  
  for( unsigned int i=0; i< Angle.size()-1; ++i ) {
    double angle =  0.5*(Angle[i]+Angle[i+1]);  // [deg]
    if(i==0 && angle_min == 0.0 ){ angle=0.0; } // first bin at z-axis
    for( unsigned int j=0; j<Energy.size()-1; j++ ){
      double energy =  0.5*(Energy[j]+Energy[j+1])*energy_scale; // [eV]
      fprintf(file, "%li  %.4e %d %.4e  %.4e %d %.4e\n", ions_tot, angle, Bin_particles[i], fact_pa_A, 
						  	      energy, Bin_particle_energy[i][j], fact_pa );
      // Angle in [degree], Energy in [eV], Bin_particles in [psudo pa], Bin_particle_energy [psudo pa], ions_tot [1]
    }
  }
  fclose(file);
  
  ions_tot = 0.; 
  std::fill( begin(Bin_particles), end(Bin_particles), 0 );
  for( auto& line : Bin_particle_energy ){
    std::fill( begin(line), end(line), 0 );
  }


#if USE_PARTICLE_ORIGIN
  const char* f_name2 = "out/ang_curr_origin.dat";
  FILE* file_origin = fopen(f_name2, "w");
  if ( !file_origin ){ 
    printf("could not open %s\n",f_name2); 
    exit(EXIT_FAILURE);
  }
  
  for( unsigned int i=0; i<Angle.size()-1; ++i) {
  double angle =  0.5*(Angle[i]+Angle[i+1]);  // [deg]
    if(i==0 && angle_min == 0.0 ){ angle=0.0; }// first bin at z-axis 
  
    fprintf(file_origin,"%.4e\n",angle); 
    for( int k=0; k<global_grid.z_dim; k++){
      for( int l=0; l<global_grid.r_dim; l++){
	fprintf(file_origin,"%d ",Bin_particle_origin[i][k][l]);
      }  
	fprintf(file_origin,"\n");
    }
  }
fclose(file_origin);
#endif

}


void out_mfp(double ecoll_dr[],double ecoll_dt[],double icoll_dr[],double icoll_dt[],double ncoll_dr[],double ncoll_dt[], std::string f_name,int n_bins )
{
  FILE    *file = fopen(f_name.c_str(), "w");
  if ( !file ){ 
    printf("could not open %s\n",f_name.c_str()); 
    exit(EXIT_FAILURE);
  }

  double fact_mfp = d_mfp*dr_0*ScaleF;	// mean free path bin [cm]
  double fact_mft = d_mft*dt_0;		// mean free time bin [s]  
  for (int i = 0; i < n_bins; ++i) { 
    fprintf(file,"%e %e %e %e %e %e %e %e \n",i*fact_mfp,i*fact_mft, ecoll_dr[i],ecoll_dt[i], icoll_dr[i],icoll_dt[i], ncoll_dr[i],ncoll_dt[i] ); 
  }
 
  fclose(file); 
}


void out_and_zero_current_to_walls( std::vector<CountParticles>& count_particles)
{
#if !NTRL_ONLY 
  double anode_current_Amp, dielectric_current_Amp, metal_current_Amp ,domain_boundary_current_Amp = 0.0;
  
  std::string f_name="out/trace_wall_current.dat";
  FILE *file_out = fopen(f_name.c_str(), "a");

  int n_aver = global_grid.layers[ELECTRONS].n_aver; 
  if(n_aver < 1){
    printf("!! n_aver=%i in out_and_zero_current_to_walls, is set to 1 for this routine\n",n_aver);
    n_aver = 1;
  }

  for (auto i = 0; i < MAX_SPECIES; ++i)
  {
    if( i == NEUTRALS ) continue;  // calculation of the current only for charged particles
    printf(" species: %s\n", (global_grid.layers[i].name).c_str());

    iprintf("n_aver=%i  I_Amp_to_n2inj=%f\n",n_aver, I_Amp_to_n2inj);
    
    // current at anode
    long int anode_charge =  global_grid.layers[i].charge * count_particles[i].anode;
    double Ianode_Amp =  (double)(anode_charge) / (n_aver * I_Amp_to_n2inj);
    anode_current_Amp += Ianode_Amp;

    // current at dielctric
    long int dielectric_charge = global_grid.layers[i].charge * count_particles[i].dielectric;
    long int dielectric_SEE_charge = global_grid.layers[i].charge * count_particles[i].dielectric_SEE;
    double Idiel_Amp = (double)(dielectric_charge-dielectric_SEE_charge) / (n_aver * I_Amp_to_n2inj);
    dielectric_current_Amp += Idiel_Amp; 

    // current at metal
    long int metal_charge = global_grid.layers[i].charge * count_particles[i].metal;
    double Imetal_Amp = (double)(metal_charge) / (n_aver * I_Amp_to_n2inj);
    metal_current_Amp += Imetal_Amp; 

    // current at domain boundary
    long int domainB_charge = global_grid.layers[i].charge * count_particles[i].domain_boundary;
    double IdomainB_Amp = (double)(domainB_charge) / (n_aver * I_Amp_to_n2inj);
    domain_boundary_current_Amp += IdomainB_Amp; 

    //debug printout
    iprintf("Total flux to walls for specie %s: charge=%i\n", (global_grid.layers[i].name).c_str(), global_grid.layers[i].charge );
    iprintf(" collected anode_charge=%li -> Ianode=%13.5e A\n", anode_charge, Ianode_Amp );
    iprintf("-> anode_current_Amp =%e\n",anode_current_Amp);
    iprintf(" dielectric_charge=%li dielectric_SEE_charge=%li -> I_diel_Amp=%13.5e A\n", dielectric_charge, dielectric_SEE_charge, Idiel_Amp);
    iprintf("-> dielectric_current_Amp=%e\n",dielectric_current_Amp);
    iprintf(" metal_charge=%li -> Imetal_Amp=%13.5e A\n", metal_charge, Imetal_Amp );
    iprintf("->metal_current_Amp=%e\n",metal_current_Amp);
    iprintf(" domainB_charge=%li -> IdomainB_Amp=%13.5e A\n", domainB_charge, IdomainB_Amp );
    iprintf("-> domain_boundary_current_Amp=%e\n",domain_boundary_current_Amp);
    
    // write in file
    fprintf(file_out,"#%s: \n",  (global_grid.layers[i].name).c_str());
    fprintf(file_out,"%li %13.5e   %li %li %13.5e   %li %13.5e   %li %13.5e \n",
		    anode_charge, Ianode_Amp, 		    
		    dielectric_charge, dielectric_SEE_charge, Idiel_Amp,
		    metal_charge, Imetal_Amp, 
		    domainB_charge, IdomainB_Amp );
  }
  fprintf(file_out, "%i %13.5e %13.5e %13.5e %13.5e \n \n", nstep, anode_current_Amp, dielectric_current_Amp, metal_current_Amp, domain_boundary_current_Amp);
  fclose(file_out);
  printf("Total currents: anode_current_Amp=%13.5eA dielectric_current_Amp=%13.5eA metal_current_Amp=%13.5eA domain_boundary_current_Amp=%13.5eA\n",
	  anode_current_Amp, dielectric_current_Amp, metal_current_Amp, domain_boundary_current_Amp );

  // zero diagnostic arrays
  for( auto& element : count_particles ){ 
    iiprintf(" zero all in element \n");
    element.anode = 0.;
    element.dielectric = 0;
    element.dielectric_SEE = 0;
    element.metal = 0.;
    element.domain_boundary =0.;
  }
#endif
}

void out_eps(vector<vector<double> > Eps)
{
 if (get_rank()==0) {
  ofstream out("eps_matrix.dat");

  for (int i = 0; i < global_grid.mesh_r_dim; ++i) {
    for (int j = 0; j < global_grid.mesh_z_dim; ++j)
    {
    out << Eps[i][j] << " ";
    }
    Eps[i].resize(global_grid.mesh_z_dim);
    out << endl;
  }
  Eps.resize(global_grid.mesh_r_dim);
 }
#if 0
  for (unsigned int r=0; r<global_grid.mesh_r_dim; ++r) {
    for (unsigned int z=0; z<global_grid.mesh_r_dim; ++z) {
    std::cout << "rank " << get_rank() << " r/z " << r << "/" << z;
    std::cout << " Eps " << Eps[r][z] << std::endl;
    }
  }
#endif
}

void out_bounds(GridLayer& layer, std::string f)
{
 if ( get_rank() == 0 ) {
  ofstream out(f+".dat");
  for (int i = 0; i < layer.r_dim; ++i){
    for (int j = 0; j < layer.z_dim; ++j){
      Cell& cell=layer.get_cell(i,j);
      out << cell.boundary_type << " ";}
      out << endl;
    }
  }
}

void out_recycle(GridLayer& layer, std::string f)
{
  if ( get_rank() == 0 ) {
    ofstream out(f+".dat");
    for (int i = 0; i < layer.r_dim; ++i){
      for (int j = 0; j < layer.z_dim; ++j){
	Cell& cell=layer.get_cell(i,j);
	out << i <<" "<<j<<" ";
	for(int rec=0;rec<cell.recycle.size();rec++){
	  out << cell.recycle[rec]<<" ";
	}
	out << endl;
      }
    }
  }
}

void out_neighbour_bounds(GridLayer& layer, std::string f)
{
 if ( get_rank() == 0 ) {
  ofstream out(f);
  for (int i = 0; i < layer.r_dim; ++i){
    for (int j = 0; j < layer.z_dim; ++j){
      Cell& cell=layer.get_cell(i,j);
      out<<i<<" "<<j<<" ";
      for(int n=0;n<8;n++){
	out<<cell.neighbours.bounds[n]<<" ";
      }
    out << endl;
    }
  out << endl;
  }
 }
}

void out_phi_bound(std::string f){
  if(!get_rank()){
    FILE *file;
    file = fopen(f.c_str(), "w");
    for(int r=0; r<global_grid.mesh_r_dim;r++){
      for(int z=0; z<global_grid.mesh_z_dim;z++){
	fprintf(file, "%f ", phi_bound[r*global_grid.mesh_z_dim+z]); 
      }
      fprintf(file, "\n");
    }
    fclose(file);
  }
}

void out_normals(GridLayer& layer, std::string f)
{
 if ( get_rank() == 0 ) {
  std::string fr= "r"+f;
  std::string fz= "z"+f;
  ofstream frout(fr);
  ofstream fzout(fz);
  for (int i = 0; i < layer.r_dim; ++i){
    for (int j = 0; j < layer.z_dim; ++j){
      Cell& cell=layer.get_cell(i,j);
      int n;
	frout << cell.normal_vector[0]<<" ";
	fzout << cell.normal_vector[1]<<" ";
      }
      frout << endl;
      fzout << endl;
    }
  }
}

void out_neighbour_normals(GridLayer& layer, std::string f)
{
 if ( get_rank() == 0 ) {
  ofstream out(f);
  for (int i = 0; i < layer.r_dim; ++i){
    for (int j = 0; j < layer.z_dim; ++j){
      Cell& cell=layer.get_cell(i,j);
      out<<i<<" "<<j<<" ";
      for(int n=0;n<8;n++){
	out<<cell.neighbours.normals[n][0]<<","<<cell.neighbours.normals[n][1]<<" ";
      }
    out << endl;
    }
  out << endl;
  }
 }
}

void type_to_color( BoundaryConditionSegment& bc_segment, ofstream& out ) {
  out << "fc rgb ";
  if ( bc_segment.type == "const_pot" ) {
    out << "\"red\" ";
    return;
  }
  if ( bc_segment.type == "dielectric_border" ) {
    out << "\"white\" ";
    return;
  }
  if ( bc_segment.type == "const_efield" ) {
    out << "\"blue\" ";
    return;
  }
  if ( bc_segment.type == "dielectric" ) {
    out << "\"white\" ";
    return;
  }
  out << "\"gray\" ";
  return;
}

void out_surf_dens(GridLayer& layer, std::string file_name)
{
  double f;
  std::vector<double> surf_dens;
  surf_dens.resize((global_grid.mesh_r_dim)*(global_grid.mesh_z_dim));

  for (int i = 0; i < layer.r_dim; ++i){
    for (int j = 0; j < layer.z_dim; ++j){
      Cell& cell=layer.get_cell(i,j);
      surf_dens[i * (layer.z_dim+1) + j]  	    += cell.surf_dens[1];
      surf_dens[i * (layer.z_dim+1) + j+1]	    += cell.surf_dens[2];
      surf_dens[(i+1) * (layer.z_dim+1) + j]	    += cell.surf_dens[0];
      surf_dens[(i+1) * (layer.z_dim+1) + j+1]	    += cell.surf_dens[3];
    }
  }

  reduce(surf_dens.data(), surf_dens.size());

 if ( get_rank() == 0 ) {
  ofstream out(file_name);
  for (int i = 0; i < layer.r_dim+1; ++i){
    for (int j = 0; j < layer.z_dim+1; ++j){
      out << surf_dens[i*(layer.z_dim+1)+j] << " ";
    }
   out << endl;
  }
 }
}

void out_pt_too_fast(GridLayer& layer, std::string file_name)
{
  double f;
  std::vector<int> too_fast;
  too_fast.resize(layer.number_of_fast_particles_vec.size());

  for (int i = 0; i < too_fast.size(); ++i){
    too_fast[i]=layer.number_of_fast_particles_vec[i];
  }
  layer.number_of_fast_particles_vec.clear();
  layer.number_of_fast_particles_vec.resize(0);

  reduce(too_fast.data(), too_fast.size());

 if ( get_rank() == 0 ) {
  ofstream out(file_name);
  for (int i = 0; i < too_fast.size(); ++i){
    out<<too_fast[i]<<std::endl;
  }
 }
}

void out_charge_dens(GridLayer& layer, std::string file_name)
{
  double f;
  std::vector<double> c_dens;
  c_dens.resize((layer.r_dim+1)*(layer.z_dim+1));
  density(layer);

  for(int i=0; i<layer.r_dim; i++){
    for(int j=0; j<layer.z_dim; j++){
      Cell& cell=layer.get_cell(i,j);
      if(i==0){
	f=-6.*qe;
	c_dens[j]  +=f*cell.dens[1]*layer.charge;
	c_dens[j+1]+=f*cell.dens[2]*layer.charge;
	//add cell.dens for grid points (1,z) too, so that the axis cell is no longer needed
	f=-qe;
	c_dens[j+(layer.z_dim+1)]  +=f*cell.dens[0]*layer.charge;
	c_dens[j+(layer.z_dim+1)+1]+=f*cell.dens[3]*layer.charge;
      } else if(i==layer.r_dim-1){
	f -= 6. * qe / (3. * layer.r_dim - 1.);  // check it out !!(is here layer.r_dim+1 needed?)
	c_dens[layer.r_dim* (layer.z_dim+1) + j]   += f * cell.dens[0]*layer.charge;
	c_dens[layer.r_dim* (layer.z_dim+1) + j+1] += f * cell.dens[3]*layer.charge;
	//add cell.dens for grid points (layer.r_dim-1,z) too, so that the axis cell is no longer needed
	f=-qe/(layer.r_dim-1);
	c_dens[i * (layer.z_dim+1) + j]	    += f * cell.dens[1]*layer.charge;  
	c_dens[i * (layer.z_dim+1) + j+1]	    += f * cell.dens[2]*layer.charge;  
      } else {
	f = -qe / i;
	c_dens[(i+1) * (layer.z_dim+1) + j]	    += f * (cell.dens[0]+cell.surf_dens[0])*layer.charge;  
	c_dens[i * (layer.z_dim+1) + j]  	    += f * (cell.dens[1]+cell.surf_dens[1])*layer.charge;  
	c_dens[i * (layer.z_dim+1) + j+1]	    += f * (cell.dens[2]+cell.surf_dens[2])*layer.charge;  
	c_dens[(i+1) * (layer.z_dim+1) + j+1]	    += f * (cell.dens[3]+cell.surf_dens[3])*layer.charge;  
      }
    }
  }

  reduce(c_dens.data(), c_dens.size());

   if ( get_rank() == 0 ) {
    ofstream out(file_name);
    for (int i = 0; i < layer.r_dim+1; ++i){
      for(int j = 0; j < layer.z_dim+1; ++j){
	out<<c_dens[i*(layer.z_dim+1)+j]<<" ";
      }
      out<<std::endl;
    }
   }
}

void type_to_facealpha( BoundaryConditionSegment& bc_segment, ofstream& out ) {
  out << "fs transparent solid ";
  if ( bc_segment.type == "const_pot" ) {
    out << "0.5";
    return;
  }
  if ( bc_segment.type == "dielectric_border" ) {
    out << "0.7";
    return;
  }
  if ( bc_segment.type == "const_efield" ) {
    out << "0.5";
    return;
  }
  if ( bc_segment.type == "dielectric" ) {
    out << "0.5";
    return;
  }
  out << "0.7";
}

/// print file with gnuplot objects coming from boundary condition file to overlay any gnuplot script
/// just add: load "bc_gnuplot.dat" 
void out_boundary_conditions_to_gnuplot( Grid& grid, BoundaryConditionSegments& bc_segments ) {

  if (get_rank() == 0) {
    ofstream out("bc_gnuplot.dat");
    for (auto& polygon : bc_segments) {
      if (polygon.point_list.size() == 0){
        continue;
      }

      out << "set object polygon from ";

      for (int j = 0; j < polygon.point_list.size(); ++j) {
        out << polygon.point_list[j].second << "," << polygon.point_list[j].first << ",1 to ";
      }

      out << polygon.point_list[0].second << "," << polygon.point_list[0].first << ",1 front ";
      type_to_color( polygon, out );
      type_to_facealpha(polygon, out );
      out << std::endl;
      }
  }
}

/// print variables used in various gnuplot scripts to avoid reading the whole slurm file
void out_variables_to_gnuplot() {

	if ( get_rank() == 0 ) {
	ofstream out("gnuvars.dat");
	out << "ScaleF=" << ScaleF << endl;
	out << "dr=" << dr_0 << endl;
	out << "dt=" << dt_0 << endl;
	double cs_help=sqrt(me_over_mi)*v_te;
	out << "cs=" << cs_help << endl;

	out << "vmax_e=" << vmax_e << endl;
	out << "dv_e=" << 2*vmax_e/NBIN<< endl;
	out << "dt_e=" << 1 << endl;

	out << "vmax_i=" << vmax_i << endl;
	out << "dv_i=" << 2*vmax_i/NBIN<< endl;
	out << "dt_i=" << dt_ion << endl;

	out << "vmax_n=" << vmax_n << endl;
	out << "dv_n=" << 2*vmax_n/NBIN<< endl;
	out << "dt_n=" << dt_ntrl << endl;
	
	out << "vmax_ni=" << vmax_ni << endl;
	out << "dv_ni=" << 2*vmax_ni/NBIN<< endl;
	out << "dt_ni=" << dt_ion << endl;

	out << "vmax_i2p=" << vmax_i2p << endl;
	out << "dv_i2p=" << 2*vmax_i2p/NBIN<< endl;
	out << "dt_i2p=" << dt_ion << endl;
	
	out << "Bcoef=" << Bcoef << endl;	
	}
}
