
#include "collisions.h"
#include "domain_decomposition.h"

void collisions(    double* engcheck,
		    Vec3d* momcheck,
		    double* coll_time,
		    double* order_time
    ) {
    
    auto begin_time = Benchmark::start();
    decompose_domain_all( global_grid );
    *order_time += Benchmark::stop(begin_time); 
    int mpi_rank = 0; 
#if USE_MPI
    MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);
#endif

    if ( mpi_rank != 0 ) {
      *engcheck = 0;
      momcheck->r = 0;
      momcheck->t = 0;
      momcheck->z = 0;

      for (int i = 0; i < 20*NBIN; ++i){
	ecoll_dr[i] = 
	  ecoll_dt[i] = 
	  ncoll_dr[i] = 
	  ncoll_dt[i] = 
	  icoll_dr[i] = 
	  icoll_dt[i] = 0;
      }

      for (int i = 0; i < NG; ++i){
	ncoll_coulomb[i] = 
	  ncoll_ion_ntrl[i] = 
	  ncoll_el_ntrl[i] = 
	  ncoll_el_ntrl_exc[i] = 
	  ncoll_ioniz[i] = 0;
      }
    }

    auto begin_coll_time = Benchmark::start();
    // collisions for Argon
#if !NTRL_ONLY
    coll_ion_ntrl(global_grid.layers[NEUTRALS], global_grid.layers[POSITIVE_IONS],
		  React_Arp_Ar, ncoll_ion_ntrl,
		  icoll_dr, icoll_dt, ncoll_dr, ncoll_dt );

    coll_coulomb(global_grid.layers[ELECTRONS], Acoll_el, momcheck,
		 engcheck, ncoll_coulomb, ecoll_dr, ecoll_dt);
    coll_el_all_fake(global_grid.layers[NEUTRALS], // molecules
		     global_grid.layers[ELECTRONS],  // electrons
		     React_Ar_el, momcheck, engcheck, ncoll_el_ntrl, 
		     ecoll_dr, ecoll_dt, ncoll_dr, ncoll_dt );
    coll_el_all_fake(global_grid.layers[NEUTRALS], // molecules
		     global_grid.layers[ELECTRONS],  // electrons
		     React_Ar_tex, momcheck, engcheck, ncoll_el_ntrl_exc,
		     ecoll_dr, ecoll_dt, ncoll_dr, ncoll_dt);

    coll_el_ntrl_ioniz(
	global_grid.layers[NEUTRALS], 
	global_grid.layers[ELECTRONS], global_grid.layers[POSITIVE_IONS], 
	React_Ar_i, momcheck, engcheck, ncoll_ioniz,
	ecoll_dr, ecoll_dt, ncoll_dr, ncoll_dt );
#endif

    reduce(ncoll_coulomb, NG);
    reduce(ncoll_ion_ntrl, NG);
    reduce(ncoll_el_ntrl, NG);
    reduce(ncoll_el_ntrl_exc, NG);
    reduce(ncoll_ioniz, NG);

    reduce( ecoll_dr, 20*NBIN);
    reduce( ecoll_dt, 20*NBIN);
    reduce( icoll_dr, 20*NBIN);
    reduce( icoll_dt, 20*NBIN);
    reduce( ncoll_dr, 20*NBIN);
    reduce( ncoll_dt, 20*NBIN);

    reduce( engcheck );
    reduce( &momcheck->r );
    reduce( &momcheck->t );
    reduce( &momcheck->z );

    *coll_time += Benchmark::stop(begin_coll_time);
}


