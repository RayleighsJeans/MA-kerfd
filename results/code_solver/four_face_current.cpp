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
												 double rP, double zP,						 // particle pos
												 double Dr, double Dz,					   // delta r/z
												 double qu,												 // charge
												 double rP0, double zP0 ) {				 // relative cell pos

	double			 current_hor1, current_vert1,
							 current_hor2, current_vert2;
	unsigned int rmax = global_grid.r_dim,
							 zmax = global_grid.z_dim;
	const char*  part;

	Cell&				 cell1 = grid_layer.get_cell(rC,zC);
	Cell&				 cell2 = grid_layer.get_cell(rC,zC);

	if ( rP0 > 0.0 ){
		if ( zP0 > 0.0 ){
			// add to current: BOTTOM_LEFT
			cell1.diagnostic_arrays.current_left += qu * Dz * (0.5+rP0P0+0.5*Dr);
			cell1.diagnostic_arrays.current_bottom += qu * Dr * (0.5+zP0P0+0.5*Dz);
			if ( (rC > 0) && (zC > 0) ) {
				Cell& cell2 = grid_layer.get_cell(rC-1,zC-1);
				cell2.diagnostic_arrays.current_top += qu * Dr * (0.5-zP0P0-0.5*Dz);
				cell2.diagnostic_arrays.current_right += qu * Dz * (0.5-rP0-0.5*Dr);
			} else if ( (rC==0) && (zC > 0) ) {
				Cell& cell2 = grid_layer.get_cell(rC,zC-1);
				cell2.diagnostic_arrays.current_bottom += qu * Dr * (0.5-zP0-0.5*Dz);
			} else if ( (rC > 0) && (zC==0) ) {
				Cell& cell2 = grid_layer.get_cell(rC-1,zC);
				cell2.diagnostic_arrays.current_left += qu * Dz * (0.5-rP0-0.5*Dr);
			}
			// diagnostic part for check of calculation
      part = "BOTTOM_LEFT"; iiprintf(">> into four face: BOTTOM_LEFT");
			current_hor1 = qu*Dz*(0.5+rP0+0.5*Dr); current_vert1 = qu*Dr*(0.5+zP0+0.5*Dz);
			current_hor2 = qu*Dz*(0.5-rP0-0.5*Dr); current_vert2 = qu*Dr*(0.5-zP0-0.5*Dz);
		} else {
			// add to current: BOTTOM_RIGHT
			cell1.diagnostic_arrays.current_right += qu * Dz * (0.5+rP0+0.5*Dr);
			cell1.diagnostic_arrays.current_bottom += qu * Dr * (0.5-zP0-0.5*Dz);
			if ( (rC < rmax-1) && (zC < zmax-1) ) {
				Cell& cell2 = grid_layer.get_cell(rC-1,zC+1);
				cell2.diagnostic_arrays.current_top += qu * Dr * (0.5+zP0+0.5*Dz);
				cell2.diagnostic_arrays.current_left += qu * Dz * (0.5-rP0-0.5*Dr);
			} else if ( (rC==0) && (zC < zmax-1) ) {
				Cell& cell2 = grid_layer.get_cell(rC,zC+1);
				cell2.diagnostic_arrays.current_bottom += qu * Dr * (0.5+zP0+0.5*Dz);
			} else if ( (rC > 0) && (zC==zmax-1) ) {
				Cell& cell2 = grid_layer.get_cell(rC-1,zC);
				cell2.diagnostic_arrays.current_right += qu * Dz * (0.5-rP0-0.5*Dr);
			}
			// diagnostic part for check of calculation
      part = "BOTTOM_RIGHT"; iiprintf(">> into four face: BOTTOM_RIGHT");
			current_hor1 = qu*Dz*(0.5+rP0+0.5*Dr); current_vert1 = qu*Dr*(0.5-zP0-0.5*Dz);
			current_hor2 = qu*Dz*(0.5-rP0-0.5*Dr); current_vert2 = qu*Dr*(0.5+zP0+0.5*Dz);
		}
	} else {								
		if ( zP0 > 0.0 ){ 
			// add to current: TOP_LEFT
			cell1.diagnostic_arrays.current_left += qu * Dz * (0.5-rP0-0.5*Dr);
			cell1.diagnostic_arrays.current_top += qu * Dr * (0.5+zP0+0.5*Dz);
			if ( (rC < rmax-1) && (zC > 0) ) {
				Cell& cell2 = grid_layer.get_cell(rC+1,zC-1);
				cell2.diagnostic_arrays.current_bottom += qu * Dr * (0.5-zP0-0.5*Dz);
				cell2.diagnostic_arrays.current_right += qu * Dz * (0.5+rP0+0.5*Dr);
			} else if ( (rC==rmax-1) && (zC > 0) ) {
				Cell& cell2 = grid_layer.get_cell(rC,zC-1);
				cell2.diagnostic_arrays.current_top += qu * Dr * (0.5-zP0-0.5*Dz);
			} else if ( (rC < rmax-1) && (zC==0) ) {
				Cell& cell2 = grid_layer.get_cell(rC+1,zC);
				cell2.diagnostic_arrays.current_left += qu * Dz * (0.5-rP0-0.5*Dr);
			}
			// diagnostic part for check of calculation
      part = "TOP_LEFT"; iiprintf(">> into four face: TOP_LEFT");
			current_hor1 = qu*Dz*(0.5-rP0-0.5*Dr); current_vert1 = qu*Dr*(0.5+zP0+0.5*Dz);
			current_hor2 = qu*Dz*(0.5+rP0+0.5*Dr); current_vert2 = qu*Dr*(0.5-zP0-0.5*Dz);
		} else {
			// add to current: TOP_RIGHT
			cell1.diagnostic_arrays.current_right += qu * Dz * (0.5-rP0-0.5*Dr);
			cell1.diagnostic_arrays.current_top += qu * Dr * (0.5-zP0-0.5*Dz);
			if ( (rC < rmax-1) && (zC < zmax-1) ) {
				Cell& cell2 = grid_layer.get_cell(rC+1,zC+1);
				cell2.diagnostic_arrays.current_bottom += qu * Dr * (0.5+zP0+0.5*Dz);
				cell2.diagnostic_arrays.current_left += qu * Dz * (0.5+rP0+0.5*Dr);
			} else if ( (rC==rmax-1) && (zC < zmax) ) {
				Cell& cell2 = grid_layer.get_cell(rC,zC+1);
				cell2.diagnostic_arrays.current_top += qu * Dr * (0.5+zP0+0.5*Dz);
			} else if ( (rC < rmax-1) && (zC==zmax-1) ) {
				Cell& cell2 = grid_layer.get_cell(rC+1,zC);
				cell2.diagnostic_arrays.current_right += qu * Dz * (0.5+rP0+0.5*Dr);
			}
			// diagnostic part for check of calculation
      part = "TOP_RIGHT"; iiprintf(">> into four face: TOP_RIGHT");
			current_hor1 = qu*Dz*(0.5-rP0-0.5*Dr); current_vert1 = qu*Dr*(0.5-zP0-0.5*Dz);
			current_hor2 = qu*Dz*(0.5+rP0+0.5*Dr); current_vert2 = qu*Dr*(0.5+zP0+0.5*Dz);
		}

		iiiprintf ( ">> four face current %s: cell indices rC=%i, zC=%i,"
							  " particle pos rP=%g, zP=%g,"
								" relavtive pos to next intersec r=%g, z=%g,"
								" second cell used rC2=%i, rZ2=%i,"
								" deltar=%g, deltaz=%g and qu=%g and finally"
								" calculated currents: cur_vert1=%g cur_hor1=%g"
								" and cur_vert2=%g cur_hor2=%g\n"
								part, rC, zC, rP, zP, r, z, cell2.r, cell2.z,
								Dr, Dz, qu, current_vert1, current_hor1,
								current_vert2, current_hor2 );

	}
}
