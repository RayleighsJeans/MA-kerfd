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
// BASE CELL CURRENT ADDITIONS
//**************************************************************

void rpzp_a2c ( rC, zC,								  // cell idices
								rP, zP,									// particle pos
								Dr, Dz,									// delta r/z
								std::string horizontal, // horizontal direction
								std::string vertical,   // vertical direction
								qu );									  // charge

void rpzm_a2c ( rC, zC,								  // cell idices
								rP, zP,									// particle pos
								Dr, Dz,									// delta r/z 
								std::string horizontal, // horizontal direction
								std::string vertical,   // vertical direction
								qu );									  // charge

void rmzp_a2c ( rC, zC,								  // cell idices
								rP, zP,									// particle pos
								Dr, Dz,									// delta r/z
								std::string horizontal, // horizontal direction
								std::string vertical,   // vertical direction
								qu );									  // charge

void rmzm_a2c ( rC, zC,								  // cell idices
								rP, zP,									// particle pos
								Dr, Dz,									// delta r/z 
								std::string horizontal, // horizontal direction
								std::string vertical,   // vertical direction
								qu );									  // charge

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

void area_density ( );

void getcellparts ( );

void cellfacecurrent ( );

//**************************************************************
// DIAGNOSTIC OUTPUT
//**************************************************************

void out_cellcurrent ( std::string dat_name );

void out_area_density ( std::string dat_name );

