

#include "collisions.h"
#include "domain_decomposition.h"

void collisions(double* engcheck, Vec3d* momcheck, double* coll_time,  double* order_time ) 
{

  auto begin_time = Benchmark::start();
  int mpi_rank = get_rank(); 

  //printf("doing other collision\n");
  decompose_domain_all( global_grid );
  *order_time += Benchmark::stop(begin_time);   


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
#if USE_TWOPLUS_IONS
	ncoll_double_ioniz[i] = 0;
	ncoll_ioniz2[i] = 0;
#endif
    }
  }
  auto begin_coll_time = Benchmark::start();
  
#if !NTRL_ONLY
  // Xe - Xep  elastic collisions (simplified version)
  coll_ion_ntrl( global_grid.layers[NEUTRALS], global_grid.layers[POSITIVE_IONS], 
		React_Xep_Xe, ncoll_ion_ntrl, icoll_dr,
		icoll_dt, ncoll_dr, ncoll_dt);

  // e - e Coulomb collisions
  coll_coulomb( global_grid.layers[ELECTRONS], Acoll_el, momcheck, engcheck, ncoll_coulomb, ecoll_dr, ecoll_dt);

  // Xe + e -> Xe + e  Total elastic collision
  coll_el_all_fake( global_grid.layers[NEUTRALS],  // molecules
		   global_grid.layers[ELECTRONS],	       // electrons
		   React_Xe_el, momcheck, engcheck, ncoll_el_ntrl,
		   ecoll_dr, ecoll_dt, ncoll_dr, ncoll_dt);

  // Xe + e -> Xe* + e  Total excitation collision
  coll_el_all_fake( global_grid.layers[NEUTRALS], // molecules
		   global_grid.layers[ELECTRONS],		       // electrons
		   React_Xe_tex, momcheck, engcheck, ncoll_el_ntrl_exc,
		   ecoll_dr, ecoll_dt, ncoll_dr, ncoll_dt);

  // Xe + e -> Xe+ + 2e   Ionisation
  coll_el_ntrl_ioniz( global_grid.layers[NEUTRALS],  
		     global_grid.layers[ELECTRONS], global_grid.layers[POSITIVE_IONS], 
		     React_Xe_i, momcheck, engcheck, ncoll_ioniz, ecoll_dr,
		     ecoll_dt, ncoll_dr, ncoll_dt);
#if USE_TWOPLUS_IONS
  // Xe + e -> Xe++ + 3e  Ionisation (2nd order)
  coll_el_ntrl_ioniz2(global_grid.layers[NEUTRALS], global_grid.layers[ELECTRONS], global_grid.layers[TWOPLUS_IONS], 
		      React_Xe_ipp, momcheck, engcheck, ncoll_double_ioniz, ecoll_dr, 
		      ecoll_dt, ncoll_dr, ncoll_dt);
 
  // Xe+ + e -> Xe++ + e   Ionisation (1st order)
  coll_el_ntrl_ioniz( global_grid.layers[POSITIVE_IONS], 
		      global_grid.layers[ELECTRONS],
		      global_grid.layers[TWOPLUS_IONS], 
		      React_Xep_ipp,
		     momcheck, engcheck, ncoll_ioniz2, ecoll_dr, ecoll_dt,
		     icoll_dr, icoll_dt);

  //TODO@julia: what about other collisions? e+i2p, i2p+i2p, i+i2p, n+i2p + excitations
#endif
#endif
 

  reduce(ncoll_coulomb, NG);
  reduce(ncoll_ion_ntrl, NG);
  reduce(ncoll_el_ntrl, NG);
  reduce(ncoll_el_ntrl_exc, NG);
  reduce(ncoll_ioniz, NG);
#if USE_TWOPLUS_IONS
  reduce(ncoll_double_ioniz, NG);
  reduce(ncoll_ioniz2, NG);
#endif

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

  *coll_time += Benchmark::stop(begin_coll_time) ;

}
