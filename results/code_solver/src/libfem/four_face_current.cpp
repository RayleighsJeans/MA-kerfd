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


void four_face_current ( unsigned int rC, unsigned int zC, // cell idices
												 double Dr, double Dz,					   // delta r/z
												 double qu,												 // charge
												 double rP0, double zP0 ) {				 // relative cell pos

  GridLayer& grid_layer = global_grid.layers[ELECTRONS];

	double			 current_hor1, current_vert1,
							 current_hor2, current_vert2;

	unsigned int rmax = global_grid.r_dim-1,
							 zmax = global_grid.z_dim-1,
               rC2, zC2;

	const char*  part;

	Cell&				 cell1 = grid_layer.get_cell(rC,zC);

  iiiprintf ( ">> before four face current: cell indices rC=%i/%i, zC=%i/%i,"
              " relavtive pos to next intersec r=%g, z=%g,"
              " deltar=%g, deltaz=%g and qu=%g\n",
              rC, rmax, zC, zmax, rP0, zP0,
              Dr, Dz, qu );

	if ( rP0 > 0.0 ){
		if ( zP0 > 0.0 ){
			// add to current: BOTTOM_LEFT
			cell1.fem.current_left += qu * Dz * (0.5+rP0+0.5*Dr);
			cell1.fem.current_bottom += qu * Dr * (0.5+zP0+0.5*Dz);
			if ( (rC > 0) && (zC > 0) ) {
				Cell& cell2 = grid_layer.get_cell(rC-1,zC-1); rC2=rC-1; zC2=zC-1;
				cell2.fem.current_top += qu * Dr * (0.5-zP0-0.5*Dz);
				cell2.fem.current_right += qu * Dz * (0.5-rP0-0.5*Dr);
			} else if ( (rC==0) && (zC > 0) ) {
				Cell& cell2 = grid_layer.get_cell(rC,zC-1); rC2=rC; zC2=zC-1;
				cell2.fem.current_bottom += qu * Dr * (0.5-zP0-0.5*Dz);
			} else if ( (rC > 0) && (zC==0) ) {
				Cell& cell2 = grid_layer.get_cell(rC-1,zC); rC2=rC-1; zC2=zC;
				cell2.fem.current_left += qu * Dz * (0.5-rP0-0.5*Dr);
			}
			// diagnostic part for check of calculation
      part = "BOTTOM_LEFT"; iiprintf(">> into four face: BOTTOM_LEFT ");
			current_hor1 = qu*Dz*(0.5+rP0+0.5*Dr); current_vert1 = qu*Dr*(0.5+zP0+0.5*Dz);
			current_hor2 = qu*Dz*(0.5-rP0-0.5*Dr); current_vert2 = qu*Dr*(0.5-zP0-0.5*Dz);
		} else {
			// add to current: BOTTOM_RIGHT
			cell1.fem.current_right += qu * Dz * (0.5+rP0+0.5*Dr);
			cell1.fem.current_bottom += qu * Dr * (0.5-zP0-0.5*Dz);
			if ( (rC > 0) && (zC < zmax) ) {
				Cell& cell2 = grid_layer.get_cell(rC-1,zC+1); rC2=rC-1; zC2=zC+1;
				cell2.fem.current_top += qu * Dr * (0.5+zP0+0.5*Dz);
				cell2.fem.current_left += qu * Dz * (0.5-rP0-0.5*Dr);
			} else if ( (rC==0) && (zC < zmax) ) {
				Cell& cell2 = grid_layer.get_cell(rC,zC+1); rC2=rC; zC2=zC+1;
				cell2.fem.current_bottom += qu * Dr * (0.5+zP0+0.5*Dz);
			} else if ( (rC > 0) && (zC==zmax) ) {
				Cell& cell2 = grid_layer.get_cell(rC-1,zC); rC2=rC-1; zC2=zC;
				cell2.fem.current_right += qu * Dz * (0.5-rP0-0.5*Dr);
			}
			// diagnostic part for check of calculation
      part = "BOTTOM_RIGHT"; iiprintf(">> into four face: BOTTOM_RIGHT ");
			current_hor1 = qu*Dz*(0.5+rP0+0.5*Dr); current_vert1 = qu*Dr*(0.5-zP0-0.5*Dz);
			current_hor2 = qu*Dz*(0.5-rP0-0.5*Dr); current_vert2 = qu*Dr*(0.5+zP0+0.5*Dz);
		}
	} else {								
		if ( zP0 > 0.0 ){ 
			// add to current: TOP_LEFT
			cell1.fem.current_left += qu * Dz * (0.5-rP0-0.5*Dr);
			cell1.fem.current_top += qu * Dr * (0.5+zP0+0.5*Dz);
			if ( (rC < rmax) && (zC > 0) ) {
				Cell& cell2 = grid_layer.get_cell(rC+1,zC-1); rC2=rC+1; zC2=zC-1;
				cell2.fem.current_bottom += qu * Dr * (0.5-zP0-0.5*Dz);
				cell2.fem.current_right += qu * Dz * (0.5+rP0+0.5*Dr);
			} else if ( (rC==rmax) && (zC > 0) ) {
				Cell& cell2 = grid_layer.get_cell(rC,zC-1); rC2=rC; zC2=zC-1;
				cell2.fem.current_top += qu * Dr * (0.5-zP0-0.5*Dz);
			} else if ( (rC < rmax) && (zC==0) ) {
				Cell& cell2 = grid_layer.get_cell(rC+1,zC); rC2=rC-1; zC2=zC;
				cell2.fem.current_left += qu * Dz * (0.5-rP0-0.5*Dr);
			}
			// diagnostic part for check of calculation
      part = "TOP_LEFT"; iiprintf(">> into four face: TOP_LEFT ");
			current_hor1 = qu*Dz*(0.5-rP0-0.5*Dr); current_vert1 = qu*Dr*(0.5+zP0+0.5*Dz);
			current_hor2 = qu*Dz*(0.5+rP0+0.5*Dr); current_vert2 = qu*Dr*(0.5-zP0-0.5*Dz);
		} else {
			// add to current: TOP_RIGHT
			cell1.fem.current_right += qu * Dz * (0.5-rP0-0.5*Dr);
			cell1.fem.current_top += qu * Dr * (0.5-zP0-0.5*Dz);
			if ( (rC < rmax) && (zC < zmax) ) {
				Cell& cell2 = grid_layer.get_cell(rC+1,zC+1); rC2=rC+1; zC2=zC+1;
				cell2.fem.current_bottom += qu * Dr * (0.5+zP0+0.5*Dz);
				cell2.fem.current_left += qu * Dz * (0.5+rP0+0.5*Dr);
			} else if ( (rC==rmax) && (zC < zmax) ) {
				Cell& cell2 = grid_layer.get_cell(rC,zC+1); rC2=rC+1; zC2=zC+1;
				cell2.fem.current_top += qu * Dr * (0.5+zP0+0.5*Dz);
			} else if ( (rC < rmax) && (zC==zmax) ) {
				Cell& cell2 = grid_layer.get_cell(rC+1,zC); rC2=rC+1; zC2=zC+1;
				cell2.fem.current_right += qu * Dz * (0.5+rP0+0.5*Dr);
			}
			// diagnostic part for check of calculation
      part = "TOP_RIGHT"; iiprintf(">> into four face: TOP_RIGHT ");
			current_hor1 = qu*Dz*(0.5-rP0-0.5*Dr); current_vert1 = qu*Dr*(0.5-zP0-0.5*Dz);
			current_hor2 = qu*Dz*(0.5+rP0+0.5*Dr); current_vert2 = qu*Dr*(0.5+zP0+0.5*Dz);
	  }
  }

		iiiprintf ( ">> four face current %s: cell indices rC=%i, zC=%i,"
								" relavtive pos to next intersec r=%g, z=%g,"
								" second cell used rC2=%i, rZ2=%i,"
								" deltar=%g, deltaz=%g and qu=%g and finally"
								" calculated currents: cur_vert1=%g cur_hor1=%g"
								" and cur_vert2=%g cur_hor2=%g\n",
								part, rC, zC, rP0, zP0, rC2, zC2,
								Dr, Dz, qu, current_vert1, current_hor1,
								current_vert2, current_hor2 );

}
