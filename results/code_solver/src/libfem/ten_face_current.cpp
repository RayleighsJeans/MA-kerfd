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
#include "fem_debug.h"

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
												double qu,												// charge
												double Dr, double Dz,   					// delta r/z
                        double rP0, double zP0 ) {        // relative particle pos in cell
	
       // rP0, zP0,            // FIRST STEP
  double	Dz1, Dr1,            // deltas for step 
					rP1, zP1,            // after first step
					Dz2, Dr2,            // SECOND STEP
          rP2, zP2,            // after second step
          Dz3, Dr3;            // THIRD STEP

  double sigr, sigz,           // quantity to distignuish
         prefix;               // move queue

  // DIVIDE INTO THREE STEPS IN HORIZONTAL,
  // VERTICAL MOVE AND THE REST
  if ( Dr > 0 ) {
    sigr = 0.5 - rP0; 
    if ( Dz > 0 ) {
      sigz = 0.5-zP0;
    } else {
      sigz = -0.5-zP0;
    }
  } else {
  sigr = -0.5 - rP0; 
    if ( Dz > 0 ) {
      sigz = 0.5-zP0;
    } else {
      sigz = -0.5-zP0;
    }
  }
 
  if        ( (rP0 + (Dr/Dz)*sigz > 0.5 ) ||
              (rP0 + (Dr/Dz)*sigz <-0.5 ) ){

      // DOING R-DIRECTION MOVE FIRST
      if ( rP0 + Dr > 0.5 ) {							    // UPWARDS MOVE
        Dr1 = 0.5-rP0;                        // delta to this move part
        rP1 =-0.5; prefix = 1.0;              // position after this move & prefix
      } else if ( rP0 + Dr <= -0.5 ) {         // DOWNWARDS MOVE
        Dr1 = -0.5-rP0;                       // delta to this move part
        rP1 = 0.5; prefix =-1.0;              // position after this move
      } else {
        iprintf(">> ten face move: no upwards/downwards move, acutally seven four face ?!\n");
        seven_face_current ( rO, zO, qu, Dr, Dz, rP0, zP0 );
        return;
      }
      // is the same for both
      Dz1 = (Dz/Dr)*Dr1;                      // corresponding move to prior
      zP1 = zP0 + Dz1;                        // other coordinate after move 
      
      // DOING R-DIRECTION SECOND
      if ( zP0 + Dz > 0.5 ) {                 // LEFT HAND MOVE
        Dz2 = 0.5-zP0-Dz1;                    // delta to this move part
        zP2 = -0.5;                           // position after this move
      } else if ( zP0 + Dz <= -0.5 ) {         // RIGHT HAND MOVE
        Dz2 = -0.5-zP0-Dz1;                   // delta to this move part
        zP2 = 0.5;                            // position after this move
      } else {
        iprintf(">> ten face move: no upwards/downwards move, acutally seven four face ?!\n");
        seven_face_current ( rO, zO, qu, Dr, Dz, rP0, zP0 );
        return;
      }
      // is the same for both
      Dr2 = (Dr/Dz)*Dz2;                      // second move delta to other coordinate
      rP2 = Dr2 - prefix*0.5;                 // second position
      Dr3 = Dr - Dr1 - Dr2;                   // final move delta 
      Dz3 = Dz - Dz1 - Dz2;                   // final move delta to other coordinate

      iiiprintf ( ">> ten face move: FIRST R, THEN Z; "
                  "sigr=%g, sigz=%g, Dr1=%g, Dz1=%g, rP1=%g, zP1=%g, "
                  "prefix=%g, Dr2=%g, Dz2=%g, rP2=%g, zP2=%g, Dr2=%g, Dz2=%g, Dr3=%g, Dz3=%g\n",
                  sigr, sigz, Dr1, Dz1, rP1, zP1, prefix, Dr2, Dz2, rP2, zP2, Dr2, Dz2, Dr3, Dz3 );

  } else if ( (zP0 + (Dz/Dr)*sigr > 0.5 ) ||
              (zP0 + (Dz/Dr)*sigr <-0.5 ) ){

      // DOING Z-DIRECTION MOVE FIRST
      if ( zP0 + Dz > 0.5 ) {							    // UPWARDS MOVE
        Dz1 = 0.5-zP0;                        // delta to this move part
        zP1 =-0.5; prefix = 1.0;              // position after this move & prefix
      } else if ( zP0 + Dz <= -0.5 ) {         // DOWNWARDS MOVE
        Dz1 = -0.5-zP0;                       // delta to this move part
        zP1 = 0.5; prefix =-1.0;              // position after this move
      } else {
        iprintf(">> ten face move: no upwards/downwards move, acutally seven four face ?!\n");
        seven_face_current ( rO, zO, qu, Dr, Dz, rP0, zP0 );
        return;
      }
      // is the same for both
      Dr1 = (Dr/Dz)*Dz1;                      // corresponding move to prior
      rP1 = rP0 + Dr1;                        // other coordinate after move 
      
      // DOING R-DIRECTION SECOND
      if ( rP0 + Dr > 0.5 ) {                 // LEFT HAND MOVE
        Dr2 = 0.5-rP0-Dr1;                    // delta to this move part
        rP2 = -0.5;                           // position after this move
      } else if ( rP0 + Dr <= -0.5 ) {         // RIGHT HAND MOVE
        Dr2 = -0.5-rP0-Dr1;                   // delta to this move part
        rP2 = 0.5;                            // position after this move
      } else {
        iprintf(">> ten face move: no left hand/ right hand move, acutally seven four face ?!\n");
        seven_face_current ( rO, zO, qu, Dr, Dz, rP0, zP0 );
        return;
      }
      // is the same for both
      Dz2 = (Dz/Dr)*Dr2;                      // second move delta to other coordinate
      zP2 = Dz2 - prefix*0.5;                 // second position
      Dz3 = Dz - Dz1 - Dz2;                   // final move delta
      Dr3 = Dr - Dr1 - Dr2;                   // final move delta to other coordinate

      iiiprintf ( ">> ten face move: FIRST Z, THEN R; "
                  "sigr=%g, sigz=%g, Dr1=%g, Dz1=%g, rP1=%g, zP1=%g, "
                  "prefix=%g, Dr2=%g, Dz2=%g, rP2=%g, zP2=%g, Dr2=%g, Dz2=%g, Dr3=%g, Dz3=%g\n",
                  sigr, sigz, Dr1, Dz1, rP1, zP1, prefix, Dr2, Dz2, rP2, zP2, Dr2, Dz2, Dr3, Dz3 );

  } else {
    waprintf( ">> ten face move: can not distignuish between step queue; discarding ...\n"
              "data: zP0+(Dz/Dr)*sigr=%g, rP0+(Dr/Dz)*sigz=%g, sigr=%g, sigz=%g\n",
              zP0+(Dz/Dr)*sigr, rP0+(Dr/Dz)*sigz, sigr, sigz );
    return;
  }

	// REAL FIRST STEP
	iiiprintf(">> ten face first step: Dr1=%g, Dz1=%g, rP0=%g, zP0=%g, rO=%d, zO=%d\n",
						Dr1, Dz1, rP0, zP0, rO, zO);
	four_face_current ( rO, zO,	 	  // cell idices
					            Dr1, Dz1,   // delta r/z
					            qu,  		 	  // charge
                      rP0, zP0 ); // relative pos in cell rO/zO

	// SECOND PART OF STEP
	iiiprintf(">> ten face second step: rO=%d, zO=%d, r2=%g, Dz2=%g"
						", rP1=%g, zP1=%g, rN1=%d, zN1=%d\n",
            rO, zO, Dr2, Dz2, rP1, zP1, (int)(rO+rP0+Dr1), (int)(zO+zP0+Dz1));
  four_face_current ( (int)(rO+rP0+Dr1), (int)(zO+zP0+Dz1),	// cell idices
                      Dr2, Dz2,                             // delta r/z
                      qu,  		 	                            // charge
                      rP1, zP1 );                           // relative pos in cell rO/zO
                      
  // THIRD PART OF STEP
  iiiprintf(">> ten face third step: Dr3=%g, Dz3=%g"
            ", rP2=%g, zP2=%g, rN2=%d, zN2=%d\n",
            Dr3, Dz3, rP2, zP2, (int)(rO+rP0+Dr2+Dr1), (int)(zO+zP0+Dz2+Dz1));
  four_face_current ( (int)(rO+rP0+Dr2+Dr1), (int)(zO+zP0+Dz2+Dz1),	  // cell idices
                      Dr3, Dz3,                                       // delta r/z
                      qu,  		 	                                      // charge
                      rP2, zP2 );                                     // relative pos in cell rO/zO

}
