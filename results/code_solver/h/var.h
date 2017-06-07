#pragma once

/**********************************************************************
                  Global variables 
 **********************************************************************/ 
#include <vector>
#include <string>

#include "extern_define.h"
#include "pic.h"

XTRN	int ne, ni, nn, nin, ni2p;

XTRN std::string geometry_file;

// permittivity matrix
XTRN std::vector<std::vector<double>>
  Eps;

XTRN    int na, 
            nXeP,
            nXe,
            n,                     
            dt_ion, ionstep, 
            dt_ntrl, ntrlstep,
            nstep,   
            maxstep,                /* end time in dt                   */
	    minstep,                /*  Begin time in dt                      */
	    step_of_biasing,
            nav_start,
            nav_time,
	    nav_dt,
	    start_flux,
	    ncoll_el, 
	    ncoll_n,
	    ncoll_n_n,
            inj_el_step,             // subcycling steps for electron injection
	    ifdst,
	    n2inj_ntrl,n2inj_ion,n2inj_el,n2inj_el_2,n2inj_el_3,n2inj_bg,
	    red_coll_z,
	    ne_pc, ni_pc, nn_pc, nin_pc,  /* RF related variables */
	    ne_tmp, ni_tmp, nn_tmp,ni2p_tmp, nin_tmp; //variables to reduce particle numbers into if USE_HDF5
	    
            
    
            
            
XTRN    double dr, dt,        /*     dr/Ldb,  dt*Omega_pe    */
        Ncell1,                 /*  Particles in first cell c      */
        n_e,                 /*  Dimensional reference electron density [cm-3]     */
        T_e,                  /*  Dimensional reference electron temperature [eV]     */ 
        Omega_ce, 	         /*  omega_ce/omega_pe        	       	    */
        qe, qi,              /*  dimensionless charges              */
        Ti_over_Te,
        dti_over_mi,		 /*  Factor fuer scaling of E_ion; */
        me_over_mi,
        me_over_mni,
        mi_over_me,
        mni_over_me,
        ma_over_me,
        mn_over_me,
	v_te,
	v_ti,
        v_tn,
        vn_thermal_reemission,
	vmax_e,
	vmax_i,
	vmax_n,
	vmax_ni,
	vmax_i2p,
        En_e,
        En_f, 
        Acoll_el, Acoll_ion,       
        n_e0, T_e0, N_SP, L_db0, Omega_pe0, dt_0, dr_0, Ua, Ua_t,Ua_sb,f_RF,
        Bcoef, 
        ScaleF,
	I_Amp_to_n2inj,          // [ number of pseudo particles / number of time steps / A ]
	I_cathode,
	I_anode,
	efield_sum,
        rec_fac, 
	asymm_fac,
	eta,
	coll_fac_ntrl_ntrl,
	pressure,
	v0_e_init,vt_e_init,
	v0_e_init_2,vt_e_init_2,
	v0_e_init_3,vt_e_init_3,
	v0_i_init,vt_i_init,
	v0_n_init,vt_n_init,vt_bg;

XTRN double collision_amplification_factor;

XTRN int atom_mass;

XTRN int BACKUP_EACH, recreate_domain_decomposition;

XTRN int print_statistics_interval;
XTRN int NNcell, NEcell;

XTRN double SEEyield;	       // yield for secondarry electron emission (SEE)

XTRN double Ampl;		//Amplification of anomalous Transport

// for HEMP 
XTRN int    Neutraliser;	// which neutraliser (TODO: neutralsier from experiment)
XTRN double Icath_Amp;          // cathode current of the neutraliser [-A]
XTRN double ne2inj_cath;	// number of injected super-electrons for the neutraliser [1]
XTRN double NeutralSource;      // propulsion source for HEMP [ real Xe0 /s ]
XTRN double SEESource_Amp;      // electron source for SEE at vessel [A]
XTRN int    SEEn2inj;		// number of injected super-electrons for source of SEE at vessel [1]

XTRN double NeutralRecyclCoef;  // to generate background pressure (maybe better with volume source)
XTRN int    old_Xenon_CX_CS; 	//flag to use wrong crossection in order to continue old runs

// DEBUG for velocities
XTRN int n_highv;


inline bool is_statistics_step() {
  return nstep%print_statistics_interval==0;
}

inline bool is_domain_recreation_step(){
  return nstep%recreate_domain_decomposition==0;
}

inline bool is_ion_step(){
  return nstep==ionstep;
}

inline bool is_ntrl_step(){
  return nstep == ntrlstep;
}
