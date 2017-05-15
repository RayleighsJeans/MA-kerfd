#include <stdio.h>
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
#include "engy.h"
#include "emission.h"
#include "mpi_wrapper.h"
#include "rng64.h"
#include "emission.h"

#include "fem_solver.h"

void getcellparts(){
#if USE_FEM_SOLVER
	
	//    |             $              |    
	//    | r+1,z       $              | r+1,z+1
	// ---|-------------$--------------|--------
	//		|             $              |
	//    |             $              |
	//    | top_left    $    top_right |
	//    |             $              |
	//    |             $              |
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//    |             $              |
	//    |             $              |
	//    |             $              |
	//    |bottom_left  $ bottom_right |
	//    |             $              |
	//    | r,z         $              | r,z+1
	// ---|-------------$--------------|-------
	//    | r-1,z       $              | r-1,z+1
	//    |             $              |

	printf(">> get old & new cell parts \n");
	unsigned int rO,zO;
	double r_old, z_old, r_new, z_new;

	for (auto& layer: global_grid.layers){ 
		if ( layer.name=="neutrals" ){ } else { 
			for ( unsigned int rN = 0; rN < layer.r_dim; ++rN ){
				for ( unsigned int zN = 0; zN < layer.z_dim; ++zN ){
					Cell& center = layer.get_cell(rN,zN);
					for ( unsigned int ptid = 0; ptid < center.size(); ++ptid ){
						Particle& pt=center.particles[ptid];
						
						r_old = pt.r_old; z_old = pt.z_old;
						r_new = pt.r; z_new = pt.z;

						//get old cell part
						//cell indices of old part
						rO = (int)r_old; zO = (int)z_old;
						if ( r_old < rO+0.5 ){ 
							if ( z_old < zO+0.5 ){
								pt.oldpart=BOTTOM_LEFT;
	    	    	} else {
								pt.oldpart=BOTTOM_RIGHT;
	    	      }
	    	    } else {								
				  		if ( z_old < zO+0.5 ){ 
								pt.oldpart=TOP_LEFT;
	    	      } else {
								pt.oldpart=TOP_RIGHT;
				  	  }
						}
						
						//get new cell part
						if ( r_new < rN+0.5 ){ 
							if ( z_new < zN+0.5 ){
								pt.newpart=BOTTOM_LEFT;
	    	    	} else {
								pt.newpart=BOTTOM_RIGHT;
	    	      }
	    	    } else {								
				  		if ( z_new < zN+0.5 ){ 
								pt.newpart=TOP_LEFT;
	    	      } else {
								pt.newpart=TOP_RIGHT;
				  	  }
						}

					} // particle ids
	 			} // zdim
			} // rdim
		} // except neutrals
	} // layers
#endif
}

void cellfacecurrent(){
#if USE_FEM_SOLVER

	getcellparts();

  //                       |r+1,z+1              + 
  //                       |                     +
  //                       |                     +
  //      bottom_right     |                     +
  //                       |                     +
  //              ptr+dz,  |     bottom_left     +
  //           ptz+dr *    |                     +
  //                 /     |                     +
  //~~~~~~~~~~~~~~~~/~~~~~~~~~~~~~~~~~~~~~       +
  //               /       |             $       +
  //        zp<---/------> |             $       +
  // ------------/-------------------------------+---------------
  //            /          |^            $       +
  //           /           ||            $       +
  //          /            ||            $       +
  //  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,+
  //        * ptr,ptz      |rp           $
  //                       |             $
  //    top_right          |    top_left $
  //                       |             $
  //                       |             $
  //                       |             $
  // r,z                   |             $
	//
	// IT IS UTTERLY IMPORTANT TO REMEMBER
	// THAT THE OLD CELL IS CORRESPONDING TO THE ORIGIN, THEREFORE
	// THE CENTER CELL IS THE DESTINATION OF FLIGHT 
	// HENCE THE OLD CELL IS NEEDED FIRST
	
	double rpO, zpO, rpN, zpN, // old and new particle pos
				 deltar, deltaz;		 // moving distance of particle

	int		 rO, zO, rN, zN,		 // old and new cell indices
				 qu;								 // charge
				 
	for (auto& layer: global_grid.layers){ 
		if ( layer.name=="neutrals" ){ } else { 
			if ( layer.name=="electrons" ){ qu=-1; } else if ( layer.name=="ions" ) { qu=1; 
			} else if ( layer.name=="nions" ) { qu=-1; } else if ( layer.name=="ions2p" ) { qu=2; } 
			for ( unsigned int rN = 0; rN < layer.r_dim; ++rN ){
				for ( unsigned int zN = 0; zN < layer.z_dim; ++zN ){
					Cell& newcell = layer.get_cell(rN,zN);
					for ( unsigned int ptid = 0; ptid < newcell.size(); ++ptid ){
						Particle& pt=newcell.particles[ptid];

						// old position
						rO = (int)pt.r_old; rpO = pt.r_old;
						zO = (int)pt.z_old; zpO = pt.z_old;
			
						// newp position
						rN = (int)pt.r; rpN = pt.r;
						zN = (int)pt.r; zpN = pt.z;

						// moving distance
						deltar=pt.r-pt.r_old;
						deltaz=pt.z-pt.z_old;
						
						// old cell
						Cell& oldcell=layer.get_cell(rO,zO);

						if ( pt.r_old < rO+0.5 ){ 
							if ( pt.z_old < zO+0.5 ){ //bottom_left
								rp=pt.r_old-rO; zp=pt.z_old-zO;
	    	    	} else {								  //bottom_right
								rp=pt.r_old-rO; zp=pt.z_old-(zO+1);
	    	      }
	    	    } else {								
				  		if ( pt.z_old < zO+0.5 ){ //top_left
								rp=pt.r_old-(rO+1); zp=pt.z_old-zO;
	    	      } else {									//top_right
								rp=pt.r_old-(rO+1); zp=pt.z_old-(zO+1);
				  	  }
						}
						
						#if 0
						  std::cout << ">> particle: oldpart=" << pt.oldpart << ", newpart=" << pt.newpart;
						  std::cout << ", r_old=" << pt.r_old << ", z_old=" << pt.z_old;
						  std::cout << ", r_new=" << pt.r << ", z_new=" << pt.z;
						  std::cout << ", rp=" << rp << ", zp=" << zp << ", deltar=" << deltar;
						  std::cout << ", deltaz=" << deltaz << ", rO=" << rO << ", zO=" << zO;
						  std::cout << ", rN=" << rN << "/" << layer.r_dim << ", zN=";
						  std::cout << zN << "/" <<  layer.z_dim << ", qu=" << qu << std::endl;
						#endif

						// ALGORITHMS FOR 4-WAY-CURRENT-ALTERATIONS
						// GOING FROM TOP_RIGHT
						if( (pt.oldpart==TOP_RIGHT) && (
								(	(rN==rO) && (zN==zO)     ) ||  //same top_right
								(	(rN==rO+1) && (zN==zO+1) ) ||  //bottom_left
								(	(rN==rO) && (zN==zO+1)   ) ||  //top_left
								(	(rN==rO+1) && (zN==zO) ) ) ) { //bottom_right

							//printf(">> transition one from top_right!\n");
							oldcell.diagnostic_arrays.current_right=qu*deltaz*(1/2-pt.r-1/2*deltar);
							oldcell.diagnostic_arrays.current_top=qu*deltar*(1/2-pt.z-1/2*deltaz);
							//4-face-transition-possible cases
							if (pt.newpart==BOTTOM_LEFT) {
								newcell.diagnostic_arrays.current_left=qu*deltaz*(1/2+pt.r+1/2*deltar);
								newcell.diagnostic_arrays.current_bottom=qu*deltar*(1/2+pt.z+1/2*deltaz);
							} else if (pt.newpart==TOP_LEFT) { 
								newcell.diagnostic_arrays.current_left=qu*deltaz*(1/2+pt.r+1/2*deltar);
								newcell.diagnostic_arrays.current_top=qu*deltar*(1/2+pt.z+1/2*deltaz);
							} else if (pt.newpart==BOTTOM_RIGHT) {
								newcell.diagnostic_arrays.current_right=qu*deltaz*(1/2+pt.r+1/2*deltar);
								newcell.diagnostic_arrays.current_bottom=qu*deltar*(1/2+pt.z+1/2*deltaz);
							} else if (pt.newpart==TOP_RIGHT) { // same cell
								newcell.diagnostic_arrays.current_right=qu*deltaz*(1/2+pt.r+1/2*deltar);
								newcell.diagnostic_arrays.current_top=qu*deltar*(1/2+pt.z+1/2*deltaz);
							}
						// GOING FROM TOP_LEFT
						} else if ( (pt.oldpart==TOP_LEFT) && (
											(	(rN==rO) && (zN==zO)     ) ||  //same top_left
										  (	(rN==rO+1) && (zN==zO-1) ) ||  //bottom_left
											(	(rN==rO) && (zN==zO-1)   ) ||  //top_left
											(	(rN==rO+1) && (zN==zO) ) ) ) { //bottom_right

							//printf(">> transition two from top_left!\n");
							oldcell.diagnostic_arrays.current_left=qu*deltaz*(1/2-pt.r-1/2*deltar);
							oldcell.diagnostic_arrays.current_top=qu*deltar*(1/2+pt.z+1/2*deltaz);
							//4-face-transition-possible cases
							if (pt.newpart==BOTTOM_LEFT) {
								newcell.diagnostic_arrays.current_left=qu*deltaz*(1/2+pt.r+1/2*deltar);
								newcell.diagnostic_arrays.current_bottom=qu*deltar*(1/2-pt.z-1/2*deltaz);
							} else if (pt.newpart==TOP_RIGHT) { 
								newcell.diagnostic_arrays.current_right=qu*deltaz*(1/2+pt.r+1/2*deltar);
								newcell.diagnostic_arrays.current_top=qu*deltar*(1/2-pt.z-1/2*deltaz);
							} else if (pt.newpart==BOTTOM_RIGHT) {
								newcell.diagnostic_arrays.current_right=qu*deltaz*(1/2+pt.r+1/2*deltar);
								newcell.diagnostic_arrays.current_bottom=qu*deltar*(1/2-pt.z-1/2*deltaz);
							} else if (pt.newpart==TOP_LEFT) { // same cell
								newcell.diagnostic_arrays.current_left=qu*deltaz*(1/2+pt.r+1/2*deltar);
								newcell.diagnostic_arrays.current_top=qu*deltar*(1/2-pt.z-1/2*deltaz);
							}
						// GOING FROM BOTTOM_LEFT
						} else if ( (pt.oldpart==BOTTOM_LEFT) && (
											(	(rN==rO) && (zN==zO)     ) ||  //same bottom_left
										  (	(rN==rO-1) && (zN==zO-1) ) ||  //top_right
											(	(rN==rO) && (zN==zO-1)   ) ||  //delta
											(	(rN==rO-1) && (zN==zO) ) ) ) { //top_left

							//printf(">> transition three from bottom_left!\n");
							oldcell.diagnostic_arrays.current_left=qu*deltaz*(1/2+pt.r+1/2*deltar);
							oldcell.diagnostic_arrays.current_bottom=qu*deltar*(1/2+pt.z+1/2*deltaz);
							//4-face-transition-possible cases
							if (pt.newpart==TOP_RIGHT) {
								newcell.diagnostic_arrays.current_right=qu*deltaz*(1/2-pt.r-1/2*deltar);
								newcell.diagnostic_arrays.current_top=qu*deltar*(1/2-pt.z-1/2*deltaz);
							} else if (pt.newpart==BOTTOM_RIGHT) { 
								newcell.diagnostic_arrays.current_right=qu*deltaz*(1/2-pt.r-1/2*deltar);
								newcell.diagnostic_arrays.current_bottom=qu*deltar*(1/2-pt.z-1/2*deltaz);
							} else if (pt.newpart==TOP_LEFT) {
								newcell.diagnostic_arrays.current_left=qu*deltaz*(1/2-pt.r-1/2*deltar);
								newcell.diagnostic_arrays.current_top=qu*deltar*(1/2-pt.z-1/2*deltaz);
							} else if (pt.newpart==BOTTOM_LEFT) { // same cell
								newcell.diagnostic_arrays.current_left=qu*deltaz*(1/2-pt.r-1/2*deltar);
								newcell.diagnostic_arrays.current_bottom=qu*deltar*(1/2-pt.z-1/2*deltaz);
							}
						// GOING FROM BOTTOM_RIGHT
						} else if ( (pt.oldpart==BOTTOM_RIGHT) && (
											(	(rN==rO) && (zN==zO)     ) ||  //same delta
										  (	(rN==rO-1) && (zN==zO+1) ) ||  //top_left
											(	(rN==rO) && (zN==zO+1)   ) ||  //bottom_left
											(	(rN==rO-1) && (zN==zO) ) ) ) { //top_right

							//printf(">> transition four from delta!\n");
							oldcell.diagnostic_arrays.current_right=qu*deltaz*(1/2+pt.r+1/2*deltar);
							oldcell.diagnostic_arrays.current_bottom=qu*deltar*(1/2-pt.z-1/2*deltaz);
							//4-face-transition-possible cases
							if (pt.newpart==TOP_LEFT) {
								newcell.diagnostic_arrays.current_left=qu*deltaz*(1/2-pt.r-1/2*deltar);
								newcell.diagnostic_arrays.current_top=qu*deltar*(1/2+pt.z+1/2*deltaz);
							} else if (pt.newpart==TOP_RIGHT) { 
								newcell.diagnostic_arrays.current_right=qu*deltaz*(1/2-pt.r-1/2*deltar);
								newcell.diagnostic_arrays.current_top=qu*deltar*(1/2+pt.z+1/2*deltaz);
							} else if (pt.newpart==BOTTOM_LEFT) {
								newcell.diagnostic_arrays.current_left=qu*deltaz*(1/2-pt.r-1/2*deltar);
								newcell.diagnostic_arrays.current_bottom=qu*deltar*(1/2+pt.z+1/2*deltaz);
							} else if (pt.newpart==BOTTOM_RIGHT) { // same cell
								newcell.diagnostic_arrays.current_right=qu*deltaz*(1/2-pt.r-1/2*deltar);
								newcell.diagnostic_arrays.current_bottom=qu*deltar*(1/2+pt.z+1/2*deltaz);
							}
						// EVERYTHING ELSE	
						} else { 
							//printf(">> there has been a more-than-4-face-crossing for current"
							//			 " or nothing happened ...\n");
						}

					} // particle ids
	 			} // zdim
			} // rdim
		} // except neutrals
	} // layers
#endif
}

