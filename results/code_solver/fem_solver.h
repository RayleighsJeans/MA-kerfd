#pragma once

#include <vector>
#include <time.h>
#include <string>
#include "var_diag.h"
#include "init_diag.h"
#include "grid_layer.h"
#include "boundary_check.h"
#include "grid.h"

//**************************************************************
// BASE CURRENT ADDITIONS
//**************************************************************
// FOUR FACE BASE CALCULATIONS
//**************************************************************

void getcellparts ( );

void four_face_current ( unsigned int rC, unsigned int zC, // cell idices
												 double rP, double zP,						 // particle pos
												 double Dr, double Dz,					   // delta r/z
												 double qu,   										 // charge
                         double rP0, double zP0 );         // relative pos in cell

void seven_face_current ( unsigned int rO, unsigned int zO, // old cell indices
													unsigned int rN, unsigned int zN, // new cell indices
													double qu,												// charge
													double Dr, double Dz,		 					// delta r/z
													double rP0, double zP0 );         // relative pos in cell

void ten_face_current (	);

void cellfacecurrent ( );

//**************************************************************
// DIAGNOSTIC OUTPUT
//**************************************************************

void out_cellcurrent ( std::string dat_name );

void out_area_density ( std::string dat_name );

//**************************************************************
// AREA CHARGE PARTS
//**************************************************************

void area_density ( );

void top_right ( GridLayer& layer, double prtr, double prtz,
		unsigned int rcenter, unsigned int zcenter,
		double Q );

void top_left ( GridLayer& layer, double prtr, double prtz,
		unsigned int rcenter, unsigned int zcenter,
		double Q );

void bottom_left ( GridLayer& layer, double prtr, double prtz,
		unsigned int rcenter, unsigned int zcenter,
		double Q );

void bottom_right ( GridLayer& layer, double prtr, double prtz,
		unsigned int rcenter, unsigned int zcenter,
		double Q );

