
#include "collisions.h"
#include "domain_decomposition.h"
#include "init_prtcls.h"

void collisions(    double* engcheck,
		    Vec3d* momcheck,
		    double* coll_time,
		    double* order_time
    ) { 

    auto begin_time = Benchmark::start();
    //printf("doing other collision\n");
    decompose_domain_all( global_grid );
    *order_time += Benchmark::stop(begin_time); 

    int mpi_rank = get_rank(); 
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
    //! collisions for Oxygen
#if USE_FEEDGAS_OXYGEN
#if !NTRL_ONLY
      //! O2+ + O2 (elastic collision)
      coll_ion_ntrl(global_grid.layers[NEUTRALS], global_grid.layers[POSITIVE_IONS], 
		    React_O2p_O2_el, ncoll_ion_ntrl,
			  icoll_dr, icoll_dt, ncoll_dr, ncoll_dt );           
      //! O- + O2 (elastic collision)
      coll_ion_ntrl(global_grid.layers[NEUTRALS], global_grid.layers[NEGATIVE_IONS], 
		    React_Omin_O2_el, ncoll_nion_ntrl,
		    icoll_dr, icoll_dt, ncoll_dr, ncoll_dt );           
      //! coulomb collision of global_grid.layers[ELECTRONS].particles
      coll_coulomb( global_grid.layers[ELECTRONS], Acoll_el,  momcheck, engcheck, ncoll_coulomb, ecoll_dr, ecoll_dt);   

      //! e + O2 (elasctic collision)
      coll_el_all_fake(global_grid.layers[NEUTRALS], 
          global_grid.layers[ELECTRONS], React_O2_el, momcheck, engcheck , 
	  ncoll_el_ntrl, ecoll_dr, ecoll_dt, ncoll_dr, ncoll_dt );
#endif
#if 0

      coll_el_all_fake(  global_grid.layers[NEUTRALS].particles,  ordcount_n, mn_over_me,  dt_ntrl,   // molecules
          global_grid.layers[ELECTRONS].particles,  ordcount_el,                     // electrons
          nr, nz,    React_O2_ex_g84, momcheck, engcheck , 
	  ncoll_el_ntrl, ecoll_dr, ecoll_dt, ncoll_dr, ncoll_dt );

      coll_el_all_fake(  global_grid.layers[NEUTRALS].particles,  ordcount_n, mn_over_me,  dt_ntrl,   // molecules
          global_grid.layers[ELECTRONS].particles,  ordcount_el,                     // electrons
          nr, nz,    React_O2_ex_g6, momcheck, engcheck , 
	  ncoll_el_ntrl, ecoll_dr, ecoll_dt, ncoll_dr, ncoll_dt );

      coll_el_all_fake(  global_grid.layers[NEUTRALS].particles,  ordcount_n, mn_over_me,  dt_ntrl,   // molecules
          global_grid.layers[ELECTRONS].particles,  ordcount_el,                     // electrons
          nr, nz,    React_O2_ex_g45, momcheck, engcheck , 
	  ncoll_el_ntrl, ecoll_dr, ecoll_dt, ncoll_dr, ncoll_dt );



      coll_el_all_fake(  global_grid.layers[NEUTRALS].particles,  ordcount_n, mn_over_me,  dt_ntrl,   // molecules
          global_grid.layers[ELECTRONS].particles,  ordcount_el,                     // electrons
          nr, nz,    React_O2_vib_1, momcheck, engcheck , 
	  ncoll_el_ntrl, ecoll_dr, ecoll_dt, ncoll_dr, ncoll_dt );

      coll_el_all_fake(  global_grid.layers[NEUTRALS].particles,  ordcount_n, mn_over_me,  dt_ntrl,   // molecules
          global_grid.layers[ELECTRONS].particles,  ordcount_el,                     // electrons
          nr, nz,    React_O2_vib_2, momcheck, engcheck , 
	  ncoll_el_ntrl, ecoll_dr, ecoll_dt, ncoll_dr, ncoll_dt );

      coll_el_all_fake(  global_grid.layers[NEUTRALS].particles,  ordcount_n, mn_over_me,  dt_ntrl,   // molecules
          global_grid.layers[ELECTRONS].particles,  ordcount_el,                     // electrons
          nr, nz,    React_O2_vib_3, momcheck, engcheck , 
	  ncoll_el_ntrl, ecoll_dr, ecoll_dt, ncoll_dr, ncoll_dt );

      coll_el_all_fake(  global_grid.layers[NEUTRALS].particles,  ordcount_n, mn_over_me,  dt_ntrl,   // molecules
          global_grid.layers[ELECTRONS].particles,  ordcount_el,                     // electrons
          nr, nz,    React_O2_vib_4, momcheck, engcheck , 
	  ncoll_el_ntrl, ecoll_dr, ecoll_dt, ncoll_dr, ncoll_dt );



      coll_el_all_fake(  global_grid.layers[NEUTRALS].particles,  ordcount_n, mn_over_me,  dt_ntrl,   // molecules
          global_grid.layers[ELECTRONS].particles,  ordcount_el,                     // electrons
          nr, nz,    React_O2_ex1D, momcheck, engcheck , 
	  ncoll_el_ntrl, ecoll_dr, ecoll_dt, ncoll_dr, ncoll_dt );

      coll_el_all_fake(  global_grid.layers[NEUTRALS].particles,  ordcount_n, mn_over_me,  dt_ntrl,   // molecules
          global_grid.layers[ELECTRONS].particles,  ordcount_el,                     // electrons
          nr, nz,    React_O2_ex1S, momcheck, engcheck , 
	  ncoll_el_ntrl, ecoll_dr, ecoll_dt, ncoll_dr, ncoll_dt );
#endif
 
#if !NTRL_ONLY
      //! e + O2 -> O2+ + 2e 
      coll_el_ntrl_ioniz( global_grid.layers[NEUTRALS], 
          global_grid.layers[ELECTRONS], global_grid.layers[POSITIVE_IONS], 
          React_O2_i,  momcheck, engcheck, 
	  ncoll_ioniz, ecoll_dr, ecoll_dt, ncoll_dr, ncoll_dt  );


      //! e + O2 -> O- + O
      coll_el_diss_attach( global_grid.layers[NEUTRALS],  
          global_grid.layers[ELECTRONS],
          global_grid.layers[NEGATIVE_IONS],
	  //global_grid.layers[NEUTRALS], 
	  ma_over_me,
          React_O2_dat,  momcheck, engcheck, 
	  ncoll_ioniz, ecoll_dr, ecoll_dt, ncoll_dr, ncoll_dt  ); 


#if 0
	

      //e + O- -> 2e + O 
      coll_el_detachment( global_grid.layers[NEGATIVE_IONS], mni_over_me, dt_ion,   // really dt_ntrl?     // Check the ion velocity!!
          global_grid.layers[ELECTRONS],
          //atoms,
          React_O_edet,  momcheck, engcheck, 
	  ncoll_ioniz, ecoll_dr, ecoll_dt, ncoll_dr, ncoll_dt) ; 


 
      //O2+ + O- ->O2 + O
      coll_ip_im_neutralization_fake(global_grid.layers[POSITIVE_IONS], mi_over_me, dt_ion,  //      // Check the ion velocity!!
          global_grid.layers[NEGATIVE_IONS], mni_over_me, dt_ion,
          React_O2_Op_Om_ntrlzFr,  momcheck, engcheck, 
	  ncoll_ioniz, ecoll_dr, ecoll_dt, ncoll_dr, ncoll_dt); 


      //e + O2 -> O + O
      coll_el_diss_recomb_fake(global_grid.layers[POSITIVE_IONS], mi_over_me, dt_ion,   // ,     // Check the ion velocity!!
          global_grid.layers[ELECTRONS],
          React_O2_dissrec,  momcheck, engcheck, 
	  ncoll_ioniz, ecoll_dr, ecoll_dt, ncoll_dr, ncoll_dt );
 
 
 
      //O- + O2 -> O + O2 + e
      coll_neutr_detachment_fake(global_grid.layers[NEUTRALS], mn_over_me, dt_ntrl,  //      // Check the ion velocity!!
          global_grid.layers[NEGATIVE_IONS], mni_over_me, dt_ion,
	  global_grid.layers[ELECTRONS], 
	  React_O2_n_detach,  momcheck, engcheck, 
	  ncoll_ioniz, ecoll_dr, ecoll_dt, ncoll_dr, ncoll_dt);  

#endif
#endif

#else

    // collisions for Argon

#if !NTRL_ONLY
    coll_ion_ntrl(global_grid.layers[NEUTRALS], global_grid.layers[POSITIVE_IONS],
		  React_Arp_Ar, ncoll_ion_ntrl,
		  icoll_dr, icoll_dt, ncoll_dr, ncoll_dt );
    coll_coulomb(global_grid.layers[ELECTRONS], Acoll_el, momcheck,
		 engcheck, ncoll_coulomb, ecoll_dr, ecoll_dt);
    coll_el_all_fake(global_grid.layers[NEUTRALS], 
		     global_grid.layers[ELECTRONS], 
		     React_Ar_el, momcheck, engcheck, ncoll_el_ntrl, ecoll_dr, ecoll_dt, ncoll_dr, ncoll_dt );
    coll_el_all_fake(global_grid.layers[NEUTRALS], 
		     global_grid.layers[ELECTRONS], 
		     React_Ar_tex, momcheck, engcheck, ncoll_el_ntrl_exc, ecoll_dr, ecoll_dt, ncoll_dr, ncoll_dt);
    coll_el_ntrl_ioniz(global_grid.layers[NEUTRALS], 
		       global_grid.layers[ELECTRONS], global_grid.layers[POSITIVE_IONS], 
		       React_Ar_i, momcheck, engcheck, ncoll_ioniz,
		       ecoll_dr, ecoll_dt, ncoll_dr, ncoll_dt );
#endif


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
    //printf(" rank %d collision time  %f\n", mpi_rank, ((double)t_coll2 - t_coll0) / CLOCKS_PER_SEC);
}

