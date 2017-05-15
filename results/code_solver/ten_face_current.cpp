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
	//              ~              $              $              
	//              ~  **************             $              
	//              ~  *           ~*             $              
	// -------~~~~~~~~~*~~~~~~~~~~~~*~~~-------------------------
	//              ~  *      *    ~*             $              
	//              ~  *   /       ~*             $
	//              ~  *  /        ~*             $              
	//          *********/***********             $              
	//          *   ~   /    *     $              $              
	//          *   ~  /     *     $              $              
	//          *   ~ *      *     $              $              
	// ---------*~~~~~~~~~~~~*------------------------------------
	//          *   ~        *                    $              
	//          *   ~        *     $              $              
	//          ****~*********     $              $              
	//              $              $              $              
	//              $              $              $              
	//              $              $              $              
	// ----------------------------------------------------------
	//              $              $              $             
	//              $              $              $             
	//              $              $              $             
	//
	// TEN FACE MOVE

void ten_face_current (	unsigned int rO, unsigned int zO, // old cell indices
												unsigned int rN, unsigned int zN, // new cell indices
												double qu,												// charge
												double Dr, double Dz,   					// delta r/z
                        double rP0, double zP0 ) {        // relative particle pos in cell
	
  double	Dz1, Dr1,            // FIRST STEP 
					rP1, zP1,            // after first step
					Dz2, Dr2,            // SECOND STEP
          rP2, zP2,            // after second step
          Dz3, Dr3;            // THIRD STEP

  // DIVIDE INTO THREE STEPS IN HORIZONTAL,
  // VERTICAL MOVE AND THE REST
	// VERTICAL
  if ( rP0 + Dr > 0.5 ) {							// UPWARDS MOVE
		Dr1 = 0.5-rP0; Dz1 = (Dr/Dz)*Dr1;
		Dz2 = Dz - Dz1; Dr2 = Dr - Dr1;
		zP1 = -0.5; zP1 = zP0 + Dz1;
	} else if ( rP0 + Dr < -0.5 ) {     // DOWNWARDS MOVE
		Dr1 = -0.5-rP0; Dz1 = (Dr/Dz)*Dr1;
		Dz2 = Dz - Dz1; Dr2 = Dr - Dr1;
		zP1 = 0.5; zP1 = zP0 + Dz1;
  }

	if ( zP0 + Dz < -0.5 ) {            // LEFT HAND MOVE
		Dz1 = -0.5-zP0; Dr1 = (Dr/Dz)*Dz1;
		Dz2 = Dz - Dz1; Dr2 = Dr - Dr1;
		zP1 = 0.5; rP1 = rP0 + Dr1;
	} else if ( zP0 + Dz > 0.5 ) {      // RIGHT HAND MOVE
		Dz1 = 0.5-zP0; Dr1 = (Dr/Dz)*Dz1;
		Dz2 = Dz - Dz1; Dr2 = Dr - Dr1;
		zP1 = -0.5; rP1 = rP0 + Dr1;
	}
	
	// REAL FIRST STEP
	iiiprintf(">> ten face first step: Dr1=%g, Dz1=%g\n"
						Dr1, Dz1);
	four_face_current ( grid_layer,// layer
	 				            rO, zO,	 	 // cell idices
					            rP0, zP0,	 // particle pos
					            Dr1, Dz1,  // delta r/z
					            qu );		 	 // charge

	// SECOND PART OF STEP
	iiiprintf(">> ten face second step: Dr2=%g, Dz2=%g"
						", rP1=%g, zP1=%g\n",	Dr2, Dz2, rP1, zP1);
	four_face_current ( grid_layer,// layer
	 				            rN, zN,		 // cell idices
					            rP1, zP1,	 // particle pos
					            Dr2, Dz2,	 // delta r/z
					            qu );			 // charge
                      
  // THIRD PART OF STEP
  iiiprintf(">> ten face third step: Dr3=%g, Dz3=%g"
            ", rP2=%g, zP2=%g\n",	Dr3, Dz3, rP2, zP2);
  four_face_current ( grid_layer,// layer
                      rN, zN,		 // cell idices
                      rP1, zP1,	 // particle pos
                      Dr2, Dz2,	 // delta r/z
                      qu );			 // charge

}
