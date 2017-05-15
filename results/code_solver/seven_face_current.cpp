#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <string>
#include <iostream>
#include <fstream>

#include "var_dim.h"
#include "var_diag.h"
#include "var.h"
#include "pic.h"
#include "arrays.h"
#include "pic.h"
#include "grid.h"
#include "emission.h"
#include "mpi_wrapper.h"

#include "fem_solver.h"
#define DEBUG_LEVEL DEBUG_INFO_3 // DEBUG_ERROR, DEBUG_INFO_1/2, DEBUG_WARNING
#include "debug_printing.h"

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
													unsigned int rN, unsigned int zN, // new cell indices
													double qu,												// charge
													double Dr, double Dz,							// delta r/z
													double rP0, double zP0 ) {        // relative pos in cell

	double	Dz1, Dr1,            // FIRST STEP 
					rP1, zP1,            // after first step
					Dz2, Dr2;            // SECOND STEP

	// if really is seven face move
	if ( ( (rP0 + Dr > 0.5) || (rP0 + Dr < -.5) ) &&
			 ( (zP0 + Dz > 0.5) || (zP0 + Dz < -.5) ) ) {
		waprintf(">> detected move over 2 inner cell boundaries; discarding seven face move ...\n");
		return;
	}

	// DECIDE WHETHER ALG IS TO BE USED
	if ( rP0 + Dr > 0.5 ) {							// UPWARDS MOVE
		Dr1 = 0.5-rP0; Dz1 = (Dr/Dz)*Dr1;
		Dz2 = Dz - Dz1; Dr2 = Dr - Dr1;
		zP1 = -0.5; zP1 = zP0 + Dz1;
	} else if ( rP0 + Dr < -0.5 ) {     // DOWNWARDS MOVE
		Dr1 = -0.5-rP0; Dz1 = (Dr/Dz)*Dr1;
		Dz2 = Dz - Dz1; Dr2 = Dr - Dr1;
		zP1 = 0.5; zP1 = zP0 + Dz1;
	} else if ( zP0 + Dz < -0.5 ) {     // LEFT HAND MOVE
		Dz1 = -0.5-zP0; Dr1 = (Dr/Dz)*Dz1;
		Dz2 = Dz - Dz1; Dr2 = Dr - Dr1;
		zP1 = 0.5; rP1 = rP0 + Dr1;
	} else if ( zP0 + Dz > 0.5 ) {      // RIGHT HAND MOVE
		Dz1 = 0.5-zP0; Dr1 = (Dr/Dz)*Dz1;
		Dz2 = Dz - Dz1; Dr2 = Dr - Dr1;
		zP1 = -0.5; rP1 = rP0 + Dr1;
	}
	
	// REAL FIRST STEP
	iiiprintf(">> seven face first step: Dr1=%g, Dz1=%g\n"
						Dr1, Dz1);
	four_face_current ( rO, zO,	 	 // cell idices
					            rP0, zP0,	 // particle pos
					            Dr1, Dz1,  // delta r/z
					            qu );		 	 // charge

	// SECOND PART OF STEP
	iiiprintf(">> seven face second step: Dr2=%g, Dz2=%g"
						", zP1=%g, rP1=%g\n",
						zP1, rP1, Dz2, Dr2);
	four_face_current ( rN, zN,		 // cell idices
					            rP1, zP1,	 // particle pos
					            Dr2, Dz2,	 // delta r/z
					            qu );			 // charge

}
