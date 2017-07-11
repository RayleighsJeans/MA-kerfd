#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string>
#include <iostream>
#include <fstream>

#include "macros.h"
#include "pic.h"
#include "arrays.h"
#include "grid.h"

#include "fem_solver.h"
#include "debug_printing.h"
//#define DEBUG_LEVEL DEBUG_ERROR

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

void getcellparts ( ) {
#if USE_FEM_SOLVER

	unsigned int		rO,zO;
	double					r_old, z_old, r_new, z_new;
	const char*			cellparts[] = {"TOP_LEFT", "TOP_RIGHT", "BOTTOM_LEFT", "BOTTOM_RIGHT"};

	for (auto& layer: global_grid.layers){ 
		if ( layer.name=="neutrals" ){ } else { 
			for ( unsigned int rN = 0; rN < layer.r_dim; ++rN ){
				for ( unsigned int zN = 0; zN < layer.z_dim; ++zN ){
					Cell& center = layer.get_cell(rN,zN);
					for ( unsigned int ptid = 0; ptid < center.size(); ++ptid ){
						Particle& pt=center.particles[ptid];

						r_old = pt.r_old; z_old = pt.z_old;
						r_new = pt.r; z_new = pt.z;
            rO = (int)r_old; zO = (int)z_old;

						//get old cell part
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

						iiiprintf ( ">> particle cell part: "
											  "ptid=%i/%i, rN=%i/%i, "
											  "zN=%i/%i, rold=%g, zold=%g"
											  " rO=%i, zO=%i"
											  " oldpart=%s, rnew=%g, znew=%g"
											  " newpart=%s\n",
											  ptid, center.size(), rN, layer.r_dim,
											  zN, layer.z_dim, r_old, z_old,
											  rO, zO, cellparts[(int)pt.oldpart], r_new, z_new,
											  cellparts[(int)pt.newpart]);

					} // particle ids
	 			} // zdim
			} // rdim
		} // except neutrals
	} // layers

#endif
}

void cellfacecurrent ( ) {
#if USE_FEM_SOLVER

	getcellparts();
	
	double			 rpO, zpO, rpN, zpN, // old and new particle pos
							 deltar, deltaz,		 // moving distance of particle
               rP0, zP0;           // relative cell pos of particle

	unsigned int rO, zO, rN, zN;		 // old and new cell indices
				double qu;								 // charge

	const char*  cellparts[] = {"TOP_LEFT", "TOP_RIGHT",
                              "BOTTOM_LEFT", "BOTTOM_RIGHT"};
				 
	for (auto& layer: global_grid.layers){ 

		// all except neutrals
		if ( layer.name=="neutrals" ){ } else { 

      // charge defintion for layer species
             if ( layer.name=="electrons" ) { qu =-1.0;
      } else if ( layer.name=="ions"			) { qu = 1.0; 
      } else if ( layer.name=="nions"			) { qu =-1.0;
      } else if ( layer.name=="ions2p"		) { qu = 2.0;
      } 

			// for layer r dim
			for ( unsigned int rN = 0; rN < layer.r_dim; ++rN ){
				for ( unsigned int zN = 0; zN < layer.z_dim; ++zN ){
					Cell& newcell = layer.get_cell(rN,zN);
					for ( unsigned int ptid = 0; ptid < newcell.size(); ++ptid ){
						Particle& pt=newcell.particles[ptid];

						// old position
						rO = (int)pt.r_old; rpO = pt.r_old;
						zO = (int)pt.z_old; zpO = pt.z_old;
			
						// new position
						rpN = pt.r;
						zpN = pt.z;

						// moving distance
						deltar=rpN-rpO;
						deltaz=zpN-zpO;

            // if no move at all
          	if ( (deltar==0.0) && (deltaz==0.0) ) {
          		waprintf(">> is not moving at all\n");
          		break;
          	} else if ( deltar == 0 ) {
          		deltar = 1e-4;
          	} else if ( deltaz == 0 ) {
          		deltaz = 1e-4;
          	}
          	
            // if more than one cell move
          	if ( sqrt(SQU(deltar)+SQU(deltaz)) >= 1.0 ) {
          		iprintf ( ">> free path=%g too big\n",
													sqrt(SQU(deltar)+SQU(deltaz)) );
							break;
          	}

						// relative position of particle in cell
          	if ( pt.oldpart==TOP_LEFT ) {
          		rP0 = rpO - (rO+1);
          		zP0 = zpO - zO;
          	} else if ( pt.oldpart==TOP_RIGHT ) {
          		rP0 = rpO - (rO+1);
          		zP0 = zpO - (zO+1);
          	} else if ( pt.oldpart==BOTTOM_LEFT ) {
          		rP0 = rpO - rO;
          		zP0 = zpO - zO;
          	} else if ( pt.oldpart==BOTTOM_RIGHT ) {
          		rP0 = rpO - rO;
          		zP0 = zpO - (zO+1);
          	}

						Cell& oldcell=layer.get_cell(rO,zO);
            
            #if 0
              // fail cause check
              // add info here at which run cacncels
              pt.oldpart = TOP_LEFT; pt.newpart = TOP_RIGHT;
              pt.r = 77.5943; pt.z = 232.519; pt.r_old = 76.6852; pt.z_old = 232.29;
              rO = 76; zO = 232; rN = 77; zN = 232; deltar = 0.909098;
              deltaz = 0.228778; qu =-1.0 ; rP0 = -0.314843; zP0 = 0.290386;
            #endif

						iiiprintf ( ">> cell current calculation: oldpart=%s, newpart=%s, r_old=%g"
												", z_old=%g, rP0=%g, zP0=%g, r_new=%g, z_new=%g, deltar=%g, deltaz=%g"
												", rO=%i, zO=%i, rN=%i/%i, zN=%i/%i, qu=%g\n",
												cellparts[(int)pt.oldpart], cellparts[(int)pt.newpart], rpO, zpO, rP0, zP0,
												rpN, zpN, deltar, deltaz, rO, zO, rN, layer.r_dim, zN, layer.z_dim, qu);

						// ALGORITHMS DECIDING BETWEEN
						// FOUR, SEVEN, AND TEN FACE CURRENT 
						// CALCULATIONS

            decide_face_move (  pt,                // particle handover
                                rO, zO,            // old cell idices
                                rN, zN,            // new cell idices
                                deltar, deltaz,		 // delta r/z
                                qu,   		 			   // charge
                                rP0, zP0 );        // relative cell pos
						
					} // particle ids
	 			} // zdim
			} // rdim
		} // except neutrals
	} // layers

#endif
}

void decide_face_move ( Particle& pt,                     // particle handover
                        unsigned int rO, unsigned int zO, // cell idices
                        unsigned int rN, unsigned int zN, // cell idices
                        double deltar, double deltaz,	    // delta r/z
                        double qu,   										  // charge
                        double rP0, double zP0 ) {        // relative cell pos
#if USE_FEM_SOLVER

      if ( 
  /***************************************************************************/
  /*********************** FOUR FACE CURRENT MOVES ***************************/ 
  /*********************** OLDPART == TOP_LEFT *******************************/ 
  /***************************************************************************/
         ( (pt.oldpart==TOP_LEFT    ) && (
         ( (pt.newpart==BOTTOM_LEFT ) && ( (rN==rO+1) && (zN==zO  ) ) ) ||
         ( (pt.newpart==BOTTOM_RIGHT) && ( (rN==rO+1) && (zN==zO-1) ) ) ||
         ( (pt.newpart==TOP_RIGHT   ) && ( (rN==rO  ) && (zN==zO-1) ) ) ||
         ( (pt.newpart==TOP_LEFT    ) && ( (rN==rO  ) && (zN==zO  ) ) ) ) )
  /***************************************************************************/
  /*********************** OLDPART == TOP_RIGHT ******************************/ 
  /***************************************************************************/
     ||  ( (pt.oldpart==TOP_RIGHT   ) && (
         ( (pt.newpart==BOTTOM_LEFT ) && ( (rN==rO+1) && (zN==zO+1) ) ) ||
         ( (pt.newpart==BOTTOM_RIGHT) && ( (rN==rO+1) && (zN==zO  ) ) ) ||
         ( (pt.newpart==TOP_RIGHT   ) && ( (rN==rO  ) && (zN==zO  ) ) ) ||
         ( (pt.newpart==TOP_LEFT    ) && ( (rN==rO  ) && (zN==zO+1) ) ) ) )
  /***************************************************************************/
  /*********************** OLDPART == BOTTOM_LEFT ****************************/ 
  /***************************************************************************/
     ||  ( (pt.oldpart==BOTTOM_LEFT ) && (
         ( (pt.newpart==BOTTOM_LEFT ) && ( (rN==rO  ) && (zN==zO  ) ) ) ||
         ( (pt.newpart==TOP_RIGHT   ) && ( (rN==rO-1) && (zN==zO-1) ) ) ||
         ( (pt.newpart==BOTTOM_RIGHT) && ( (rN==rO  ) && (zN==zO-1) ) ) ||
         ( (pt.newpart==TOP_LEFT    ) && ( (rN==rO-1) && (zN==zO  ) ) ) ) )
  /***************************************************************************/
  /*********************** OLDPART == BOTTOM_RIGHT ***************************/ 
  /***************************************************************************/
     ||  ( (pt.oldpart==BOTTOM_RIGHT) && (
         ( (pt.newpart==BOTTOM_LEFT ) && ( (rN==rO  ) && (zN==zO+1) ) ) ||
         ( (pt.newpart==TOP_RIGHT   ) && ( (rN==rO-1) && (zN==zO  ) ) ) ||
         ( (pt.newpart==BOTTOM_RIGHT) && ( (rN==rO  ) && (zN==zO  ) ) ) ||
         ( (pt.newpart==TOP_LEFT    ) && ( (rN==rO-1) && (zN==zO+1) ) ) ) )
       ) { four_face_current ( rO, zO, deltar, deltaz, qu, rP0, zP0 );
  } else if ( 
  /***************************************************************************/
  /*********************** SEVEN FACE CURRENT MOVES **************************/ 
  /************************* OLDPART == TOP_LEFT *****************************/ 
  /***************************************************************************/
         ( (pt.oldpart==TOP_LEFT    ) && (
         ( (pt.newpart==BOTTOM_LEFT ) && ( ( (rN==rO  ) && (zN==zO  ) )   ||
                                           ( (rN==rO+1) && (zN==zO+1) )   ||
                                           ( (rN==rO+1) && (zN==zO-1) ) ) )
      || ( (pt.newpart==TOP_RIGHT   ) && ( ( (rN==rO  ) && (zN==zO  ) )   ||
                                           ( (rN==rO-1) && (zN==zO-1) )   ||
                                           ( (rN==rO+1) && (zN==zO-1) ) ) )
      || ( (pt.newpart==BOTTOM_RIGHT) && ( ( (rN==rO  ) && (zN==zO-1) )   ||
                                           ( (rN==rO+1) && (zN==zO  ) ) ) ) 
      || ( (pt.newpart==TOP_LEFT    ) && ( ( (rN==rO+1) && (zN==zO  ) )   ||
                                           ( (rN==rO  ) && (zN==zO-1) )   ||
                                           ( (rN==rO  ) && (zN==zO+1) )   ||
                                           ( (rN==rO-1) && (zN==zO  ) ) ) ) ) )
  /***************************************************************************/
  /************************* OLDPART == TOP_RIGHT ****************************/ 
  /***************************************************************************/
     ||  ( (pt.oldpart==TOP_RIGHT   ) && (
         ( (pt.newpart==BOTTOM_LEFT ) && ( ( (rN==rO  ) && (zN==zO+1) )   ||
                                           ( (rN==rO+1) && (zN==zO  ) ) ) )
      || ( (pt.newpart==TOP_RIGHT   ) && ( ( (rN==rO  ) && (zN==zO-1) )   ||
                                           ( (rN==rO-1) && (zN==zO  ) )   ||
                                           ( (rN==rO  ) && (zN==zO+1) )   ||
                                           ( (rN==rO+1) && (zN==zO  ) ) ) )
      || ( (pt.newpart==BOTTOM_RIGHT) && ( ( (rN==rO  ) && (zN==zO  ) )   ||
                                           ( (rN==rO+1) && (zN==zO-1) )   ||
                                           ( (rN==rO+1) && (zN==zO+1) ) ) )
      || ( (pt.newpart==TOP_LEFT    ) && ( ( (rN==rO  ) && (zN==zO  ) )   ||
                                           ( (rN==rO-1) && (zN==zO+1) )   ||
                                           ( (rN==rO+1) && (zN==zO+1) ) ) ) ) )
  /***************************************************************************/
  /************************* OLDPART == BOTTOM_LEFT **************************/ 
  /***************************************************************************/
     ||  ( (pt.oldpart==BOTTOM_LEFT ) && (
         ( (pt.newpart==BOTTOM_LEFT ) && ( ( (rN==rO-1) && (zN==zO  ) )   ||
                                           ( (rN==rO  ) && (zN==zO-1) )   ||
                                           ( (rN==rO  ) && (zN==zO+1) )   ||
                                           ( (rN==rO+1) && (zN==zO  ) ) ) )
      || ( (pt.newpart==TOP_RIGHT   ) && ( ( (rN==rO-1) && (zN==zO  ) )   ||
                                           ( (rN==rO  ) && (zN==zO-1) ) ) )
      || ( (pt.newpart==BOTTOM_RIGHT) && ( ( (rN==rO  ) && (zN==zO  ) )   ||
                                           ( (rN==rO-1) && (zN==zO-1) )   ||
                                           ( (rN==rO+1) && (zN==zO-1) ) ) ) 
      || ( (pt.newpart==TOP_LEFT    ) && ( ( (rN==rO  ) && (zN==zO  ) )   ||
                                           ( (rN==rO-1) && (zN==zO+1) )   ||
                                           ( (rN==rO-1) && (zN==zO-1) ) ) ) ) )
  /***************************************************************************/
  /************************* OLDPART == BOTTOM_RIGHT **************************/ 
  /***************************************************************************/
     ||  ( (pt.oldpart==BOTTOM_RIGHT) && (
         ( (pt.newpart==BOTTOM_LEFT ) && ( ( (rN==rO  ) && (zN==zO+1) )   ||
                                           ( (rN==rO-1) && (zN==zO+1) )   ||
                                           ( (rN==rO+1) && (zN==zO+1) ) ) )
      || ( (pt.newpart==TOP_RIGHT   ) && ( ( (rN==rO-1) && (zN==zO-1) )   ||
                                           ( (rN==rO-1) && (zN==zO+1) )   ||
                                           ( (rN==rO  ) && (zN==zO  ) ) ) )
      || ( (pt.newpart==BOTTOM_RIGHT) && ( ( (rN==rO  ) && (zN==zO+1) )   ||
                                           ( (rN==rO  ) && (zN==zO-1) )   ||
                                           ( (rN==rO+1) && (zN==zO  ) )   ||
                                           ( (rN==rO-1) && (zN==zO  ) ) ) )
      || ( (pt.newpart==TOP_LEFT    ) && ( ( (rN==rO  ) && (zN==zO+1) )   ||
                                           ( (rN==rO-1) && (zN==zO  ) ) ) ) ) )
       ) { seven_face_current ( rO, zO, qu, deltar, deltaz, rP0, zP0 );
  /***************************************************************************/
  /****************** EVERYTHING ELSE IS TEN FACE CURRENT MOVE ***************/
  /***************************************************************************/
  } else { ten_face_current ( rO, zO, qu, deltar, deltaz, rP0, zP0 );         }

#endif
}
