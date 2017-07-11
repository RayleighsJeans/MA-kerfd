#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <string>
#include <iostream>
#include <fstream>

#include "grid.h"

#include "fem_solver.h"
#include "debug_printing.h"
//#define DEBUG_LEVEL DEBUG_ERROR

	//              $              $              $              
	//              $              $              $              
	//              $              $              $              
	//              $              $              $              
	//              $              $              $              
	// ----------------------------------------------------------
	//              $              $              $              
	//              $       *************         $              
	//              $       *      ~    *         $              
	//          **************     ~    *         $              
	//          *   ~       **     *    *         $              
	//          *   ~       **     ~    *         $              
	//          *   ~ *     **     ~    *         $              
	// ---------*~~~~~~~~~~~**~~~~~~~~~~*------------------------
	//          *   ~       *************         $              
	//          *   ~        *     $              $              
	//          **************     $              $              
	//              $              $              $              
	//              $              $              $              
	//              $              $              $              
	// ----------------------------------------------------------
	//              $              $              $             
	//              $              $              $             
	//              $              $              $             
	//
	// SEVEN FACE MOVE
	// note that there are seven cell faces touched by the unit square 
	// charge area: 3 for each position, end and start, and one both share
	// this is a necessity! 
	// there are four possible outcomes, differing by which way the charge moves from
	// its old position: upwards, downwards, left and right
	// for seven face-movements it is only needed to search for
	// rP +/- deltar >=/<= 0.5/-0.5 
	// OR
	// zP +/- deltaz >=/<= 0.5/-0.5 
	// so the charge moves over one inner cell part boundary
	// according to this there are 4 different calculation methods
	// where the 2nd part/2nd step of the move is almost the same for 
	// all of them

void seven_face_current ( unsigned int rO, unsigned int zO, // old cell indices
													double qu,												// charge
													double Dr, double Dz,							// delta r/z
													double rP0, double zP0 ) {        // relative pos in cell
#if USE_FEM_SOLVER

	double	Dz1, Dr1,            // FIRST STEP 
					rP1, zP1,            // after first step
					Dz2, Dr2;            // SECOND STEP

	// if really is seven face move
	if ( ( (rP0 + Dr > 0.5) || (rP0 + Dr < -.5) ) &&
			 ( (zP0 + Dz > 0.5) || (zP0 + Dz < -.5) ) ) {
		waprintf(">> detected move over 2 inner cell boundaries;"
             " discarding seven face move for ten face ...\n");
    ten_face_current ( rO, zO, qu, Dr, Dz, rP0, zP0 );
		return;
	}

	// DECIDE WHICH ALG IS TO BE USED
	if ( rP0 + Dr > 0.5 ) {							// UPWARDS MOVE
		Dr1 = 0.5-rP0; Dz1 = (Dz/Dr)*Dr1;
		Dz2 = Dz - Dz1; Dr2 = Dr - Dr1;
		rP1 = -0.5; zP1 = zP0 + Dz1;
	} else if ( rP0 + Dr < -0.5 ) {     // DOWNWARDS MOVE
		Dr1 = -0.5-rP0; Dz1 = (Dz/Dr)*Dr1;
		Dz2 = Dz - Dz1; Dr2 = Dr - Dr1;
		rP1 = 0.5; zP1 = zP0 + Dz1;
	} else if ( zP0 + Dz < -0.5 ) {     // LEFT HAND MOVE
		Dz1 = -0.5-zP0; Dr1 = (Dr/Dz)*Dz1;
		Dz2 = Dz - Dz1; Dr2 = Dr - Dr1;
		zP1 = 0.5; rP1 = rP0 + Dr1;
	} else if ( zP0 + Dz > 0.5 ) {      // RIGHT HAND MOVE
		Dz1 = 0.5-zP0; Dr1 = (Dr/Dz)*Dz1;
		Dz2 = Dz - Dz1; Dr2 = Dr - Dr1;
		zP1 = -0.5; rP1 = rP0 + Dr1;
	} else {
    waprintf(">> neither way move, check for four face move ...\n");
    four_face_current ( rO, zO, qu, Dr, Dz, rP0, zP0 );
    return;
  }
	
	// REAL FIRST STEP
	iiiprintf(">> seven face first step: Dr1=%g, Dz1=%g\n",
						Dr1, Dz1);
	four_face_current ( rO, zO,	 	 // cell idices
					            Dr1, Dz1,  // delta r/z
					            qu, 		 	 // charge
                      rP0, zP0 );// particle pos

	// SECOND PART OF STEP
	iiiprintf(">> seven face second step: Dr2=%g, Dz2=%g"
						", zP1=%g, rP1=%g, zN=%d, rN=%d\n",
						Dr2, Dz2, zP1, rP1, (int)(zO+zP0+Dz), (int)(rO+rP0+Dr));
	four_face_current ( (int)(rO+rP0+Dr), (int)(zO+zP0+Dz), // cell idices
					            Dr2, Dz2,	                          // delta r/z
					            qu, 			                          // charge
                      rP1, zP1 );                         // particle pos

#endif
}
