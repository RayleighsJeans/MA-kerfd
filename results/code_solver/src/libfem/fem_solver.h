#pragma once

#include <vector>
#include <time.h>
#include <string>
#include "grid_layer.h"
#include "grid.h"

//**************************************************************
// BASE CURRENT ADDITIONS
//**************************************************************
// MAIN PARTS
//**************************************************************

void getcellparts ( );

void cellfacecurrent ( );

void decide_face_move ( Particle& pt,                     // particle handover
                        unsigned int rO, unsigned int zO, // cell idices
                        unsigned int rN, unsigned int zN, // cell idices
                        double deltar, double deltaz,	    // delta r/z
                        double qu,   										  // charge
                        double rP0, double zP0 );         // relative cell pos

//**************************************************************
// CURRENT FACE CALCULATIONS 
//**************************************************************

void four_face_current ( unsigned int rC, unsigned int zC, // cell idices
												 double Dr, double Dz,					   // delta r/z
												 double qu,   										 // charge
                         double rP0, double zP0 );         // relative pos in cell

void seven_face_current ( unsigned int rO, unsigned int zO, // old cell indices
													double qu,												// charge
													double Dr, double Dz,		 					// delta r/z
													double rP0, double zP0 );         // relative pos in cell

void ten_face_current (	unsigned int rO, unsigned int zO, // old cell indices
                        double qu,												// charge
                        double Dr, double Dz,   					// delta r/z
                        double rP0, double zP0 );         // relative particle pos in cell

//**************************************************************
// ELECTRIC FIELD CALCULATIONS
//**************************************************************

void efield_matrix2fem ( Field E_grid[] );

void   efield_constpot ( unsigned int r, unsigned int z,
                         std::vector<int>& epsvec, double Ua_SB,
                         Cell& cell, BoundarySegment& polygon,
                         std::vector<bool> &queue );

void efield_constfield ( unsigned int r, unsigned int z,
                         std::vector<int>& epsvec,
                         Cell& cell, BoundarySegment& polygon,
                         std::vector<bool> &queue );

void fem_eps2cell ( );

std::vector<int>
         eps_grid2face ( unsigned int r, unsigned int z,
                         unsigned int rmax, unsigned int zmax );

void efield_fem ( );

void store_old_field ( );

void reset_fem ( ) ;

//**************************************************************
// DIAGNOSTIC OUTPUT
//**************************************************************

void output_fem ( Field E_grid[] );

void out_cellcurrent ( std::string dat_name );

void out_area_density ( std::string dat_name );

void out_efield ( std::string dat_name1, std::string dat_name2,
                  std::string dat_name3, Field E_grid[],
                  std::string dat_name4 );
 
void out_vector_femfield ( std::string dat_name );

void out_vector_stdfield ( Field E_grid[], std::string dat_name );

std::vector<std::vector<double>> out_std_efield ( Field E_grid[], std::string dat_name2 );

//**************************************************************
// AREA CHARGE PARTS
//**************************************************************

void area_density ( );

void top_right ( double prtr, double prtz,
		unsigned int rcenter, unsigned int zcenter,
		double Q, unsigned int n_aver  );

void top_left ( double prtr, double prtz,
		unsigned int rcenter, unsigned int zcenter,
		double Q, unsigned int n_aver  );

void bottom_left ( double prtr, double prtz,
		unsigned int rcenter, unsigned int zcenter,
		double Q, unsigned int n_aver  );

void bottom_right ( double prtr, double prtz,
		unsigned int rcenter, unsigned int zcenter,
		double Q, unsigned int n_aver  );
