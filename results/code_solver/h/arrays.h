#pragma once

/************************************************************************************

			      Definition der Datenfelder


 ************************************************************************************/
#include "extern_define.h"
#include "pic.h"
#include <vector>

XTRN Field       *E_grid, 
                 *E_grid_old, 
		 *E_ion,
                 *B_el, *B_ion;
                
                 
XTRN double      *phi, *phi_av,
		 *phi_bound,
		 *phi_old,
                 *dens_el,
                 *dens_ion, 
                 *dens_ion2p, 
                 *dens_nion;  

//collision diagnostics 
XTRN int      	 *ncoll_coulomb, *ncoll_ioniz, *ncoll_ntrl_ntrl, *ncoll_ion_ntrl, *ncoll_nion_ntrl, *ncoll_el_ntrl, *ncoll_el_ntrl_exc;
XTRN int 	 *ncoll_double_ioniz;   // for:  n0 + e- -> i++ + 2e- 
XTRN int  	 *ncoll_ioniz2; 	// for:  i+ + e- -> i++ + e-
XTRN int      	 *ncoll_dat, *ncoll_ntrlz, *ncoll_recomb, *ncoll_el_det, *ncoll_ntrl_detach;

XTRN double    	 *ecoll_dr, *ecoll_dt, *icoll_dr, *icoll_dt, *nicoll_dr, *nicoll_dt, *ncoll_dr, *ncoll_dt;

//for domain decomposition
XTRN int         *ordcount_el, *ordind_el, *ordcount_n, *ordind_n, 
		 *ordcount_i, *ordind_i, *ordcount_ni, *ordind_ni, *ordcount_i2p, *ordind_i2p, 
		 *ordcount_a; 

//XTRN char         *gone_n;                               
                
XTRN double      *surf_densR, *surf_densZ;  

XTRN int         *recycleR, *recycleZ;                
XTRN int         *recycleR1, *recycleR2, *recycleA, *recycleE, *recycleDFront, *recycleM;
                 
XTRN Moments     *mom_ion, *mom_ion2p,
		 *mom_nion,
                 *mom_el, 
                 *mom_n;       
                 
XTRN int         *Pion;
