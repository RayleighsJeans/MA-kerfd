#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <string>
#include <iostream>
#include <fstream>

#include "var_dim.h"
#include "density.h"
#include "var_diag.h"
#include "var.h"
#include "pic.h"
#include "arrays.h"
#include "pic.h"
#include "grid.h"
#include "init_diag.h"
#include "macros.h"
#include "output.h"
#include "engy.h"
#include "emission.h"
#include "mpi_wrapper.h"
#include "rng64.h"
#include "emission.h"

#include "fem_solver.h"

void alpha( GridLayer& layer, double prtr, double prtz,
		unsigned int rcenter, unsigned int zcenter,
		double Q ){
#if USE_FEM_SOLVER
	// printf("using TOP_LEFT; \n");
	// using algorithm delta
	double eta=(zcenter+2)-(prtz-0.5); 
	double ksi=prtr+0.5-rcenter-1;       

	// get relevant neighbouring cells 
	// CENTER, NORTH, WEST, NORTHWEST
	// calculating area charge for four cells successively
	
	// CENTER
	Cell& middle = layer.get_cell(rcenter,zcenter);
	double fi = 1./(Ncell1*(rcenter+0.5)*layer.n_aver);
	middle.diagnostic_arrays.area_weighted_charge+=Q*fi*((1.0-ksi)*(1.0-eta));
	// printf("assigned to CENTER: %g, area=%g\n", Q*fi*((1.0-eta)*(1.0-ksi)), (1.0-eta)*(1.0-ksi));
	// printf(", charge %g, fi=%g, eta=%g, ksi=%g \n", Q, fi, eta, ksi);

	// NORTHWEST
	if ( (rcenter < layer.r_dim-2) && (zcenter > 0) ) {
		fi = 1./(Ncell1*(rcenter+1+0.5)*layer.n_aver);
		// printf(", to NW: %g, area=%g\n", Q*fi*(eta*ksi), (eta*ksi));
		// printf(", rNW=%i zNW=%i\n", rcenter+1,zcenter-1);
		int rNW=rcenter+1; int zNW=zcenter-1; Cell& northwest = layer.get_cell(rNW,zNW);
		northwest.diagnostic_arrays.area_weighted_charge+=Q*fi*(eta*ksi);
	}

	// WEST
	if ( zcenter > 0 ) {
		fi = 1./(Ncell1*(rcenter+0.5)*layer.n_aver);
		// printf(", to W: %g, area=%g\n", Q*fi*(eta*(1.0-ksi)), (eta*(1.0-ksi)));
		// printf(", rW=%i zW=%i\n", rcenter,zcenter-1);
		int rW=rcenter; int zW=zcenter-1;	Cell& west = layer.get_cell(rW,zW);
		west.diagnostic_arrays.area_weighted_charge+=Q*fi*(eta*(1.0-ksi));
	}

	// NORTH
	if ( rcenter < layer.r_dim-2 ) {
		fi = 1./(Ncell1*(rcenter+1+0.5)*layer.n_aver);
		// printf(", to N: %g, area=%g\n", Q*fi*(ksi*(1.0-eta)), ((1.0-eta)*ksi));
		// printf(", rN=%i zN=%i\n", rcenter+1,zcenter);
		int rN=rcenter+1; int zN=zcenter;	Cell& north = layer.get_cell(rN,zN);
		north.diagnostic_arrays.area_weighted_charge+=Q*fi*(ksi*(1.0-eta));
	}
#endif
}

void beta( GridLayer& layer, double prtr, double prtz,
		unsigned int rcenter, unsigned int zcenter,
		double Q ){
#if USE_FEM_SOLVER
	// printf("using TOP_RIGHT; \n");
	// using algorithm gamma
	double eta=prtz+0.5-zcenter-1;
	double ksi=prtr+0.5-rcenter-1;

	// get relevant neighbouring cells 
	// NORTH, WEST, NORTHWEST
	// calculating area charge for four cells successively
	
	// CENTER
	Cell& middle = layer.get_cell(rcenter,zcenter);
	double fi = 1./(Ncell1*(rcenter+0.5)*layer.n_aver);
	middle.diagnostic_arrays.area_weighted_charge+=Q*fi*((1.0-ksi)*(1.0-eta));
	// printf("assigned to CENTER: %g, area=%g\n", Q*fi*((1.0-eta)*(1.0-ksi)), (1.0-eta)*(1.0-ksi));
	// printf(", charge %g, fi=%g, eta=%g, ksi=%g \n", Q, fi, eta, ksi);
	
	// NORTHEAST
	if ( (rcenter < layer.r_dim-2) && (zcenter < layer.z_dim-2) ) {
		fi = 1./(Ncell1*(rcenter+1+0.5)*layer.n_aver);
		// printf(", to NE: %g, area=%g\n", Q*fi*(eta*ksi), (eta*ksi));
		// printf(", rNE=%i zNE=%i\n", rcenter+1,zcenter+1);
		int rNE=rcenter+1; int zNE=zcenter+1; Cell& northeast = layer.get_cell(rNE,zNE);
		northeast.diagnostic_arrays.area_weighted_charge+=Q*fi*(eta*ksi);
	}
	
	// EAST
	if ( zcenter < layer.z_dim-2 ) {
		fi = 1./(Ncell1*(rcenter+0.5)*layer.n_aver);
		// printf(", to E: %g, area=%g\n", Q*fi*(eta*(1.0-ksi)), (eta*(1.0-ksi)));
		// printf(", rE=%i zE=%i\n", rcenter,zcenter+1);
		int rE=rcenter; int zE=zcenter+1;	Cell& east = layer.get_cell(rE,zE);
		east.diagnostic_arrays.area_weighted_charge+=Q*fi*(eta*(1.0-ksi));
	}

	// NORTH
	if ( rcenter < layer.r_dim-2 ) {
		fi = 1./(Ncell1*(rcenter+0.5)*layer.n_aver);
		// printf(", to N: %g, area=%g\n", Q*fi*((1.0-eta)*ksi), ((1.0-eta)*ksi));
		// printf(", rN=%i zN=%i\n", rcenter+1,zcenter);
		int rN=rcenter+1; int zN=zcenter;	Cell& north = layer.get_cell(rN,zN);
		north.diagnostic_arrays.area_weighted_charge+=Q*fi*(ksi*(1.0-eta));
	}
#endif
}

void gamma( GridLayer& layer, double prtr, double prtz,
		unsigned int rcenter, unsigned int zcenter,
		double Q ){
#if USE_FEM_SOLVER
	// printf("using BOTTOM_LEFT; \n");
	// using algorithm beta
	double ksi=(rcenter-1+1)-(prtr-0.5);
	double eta=zcenter-1+1-(prtz-0.5);

	// get relevant neighbouring cells 
	// WEST, SOUTH, SOUTH-WEST
	// calculating area charge for four cells successively
	
	// CENTER
	Cell& middle = layer.get_cell(rcenter,zcenter);
	double fi = 1./(Ncell1*(rcenter+0.5)*layer.n_aver);
	middle.diagnostic_arrays.area_weighted_charge+=Q*fi*((1.0-ksi)*(1.0-eta));
	// printf("assigned to CENTER: %g, area=%g\n", Q*fi*((1.0-eta)*(1.0-ksi)), (1.0-eta)*(1.0-ksi));
	// printf(", charge %g, fi=%g, eta=%g, ksi=%g \n", Q, fi, eta, ksi);

	// SOUTH-WEST
	if ( (rcenter > 0) && (zcenter > 0) ) {
		fi = 1./(Ncell1*(rcenter-1+0.5)*layer.n_aver);
		// printf(", to SW: %g, area=%g\n", Q*fi*(eta*ksi), (eta*ksi));
		// printf(", rS=%i zS=%i\n", rcenter-1,zcenter-1);
		int rSW=rcenter-1; int zSW=zcenter-1; Cell& southwest = layer.get_cell(rSW,zSW);
		southwest.diagnostic_arrays.area_weighted_charge+=Q*fi*(eta*ksi);
	}
	
	// WEST
	if ( zcenter > 0 ) {
		fi = 1./(Ncell1*(rcenter+0.5)*layer.n_aver);
		// printf(", to W: %g, area=%g\n", Q*fi*(eta*(1.0-ksi)), (eta*(1.0-ksi)));
		// printf(", rW=%i zW=%i\n", rcenter,zcenter-1);
		int rW=rcenter; int zW=zcenter-1; Cell& west = layer.get_cell(rW,zW);
		west.diagnostic_arrays.area_weighted_charge+=Q*fi*(eta*(1.0-ksi));
	}
	
	// SOUTH
	if ( rcenter > 0 ) {
		fi = 1./(Ncell1*(rcenter-1+0.5)*layer.n_aver);
		// printf(", to S: %g, area=%g\n", Q*fi*((1.0-eta)*ksi), ((1.0-eta)*ksi));
		// printf(", rS=%i zS=%i\n", rcenter-1,zcenter);
		int rS=rcenter-1; int zS=zcenter; Cell& south = layer.get_cell(rS,zS);
		south.diagnostic_arrays.area_weighted_charge+=Q*fi*(ksi*(1.0-eta));
	}
#endif
}

void delta( GridLayer& layer, double prtr, double prtz,
		unsigned int rcenter, unsigned int zcenter,
		double Q ){
#if USE_FEM_SOLVER
  // printf("using BOTTOM_RIGHT; \n");
	// using algorithm alpha
	double eta=prtz+0.5-zcenter-1;
	double ksi=(rcenter-1+1)-(prtr-0.5);

 	// get relevant neighbouring cells 
	// calculating area charge for four cells successively
	// EAST, SOUTH, SOUTH-EAST
	
	// CENTER
	Cell& middle = layer.get_cell(rcenter,zcenter);
	double fi = 1./(Ncell1*(rcenter+0.5)*layer.n_aver);
	middle.diagnostic_arrays.area_weighted_charge+=Q*fi*((1.0-ksi)*(1.0-eta));
	// printf("assigned to CENTER: %g, area=%g\n", Q*fi*((1.0-eta)*(1.0-ksi)), (1.0-eta)*(1.0-ksi));
	// printf(", charge %g, fi=%g, eta=%g, ksi=%g \n", Q, fi, eta, ksi);

	// SOUTH-EAST
	if ( (rcenter > 0) && (zcenter < layer.z_dim-2) ) {
		fi = 1./(Ncell1*(rcenter-1+0.5)*layer.n_aver);
		// printf(", to SE: %g, area=%g\n", Q*fi*(eta*ksi), (eta*ksi));
		// printf(", rSE=%i zSE=%i\n", rcenter-1,zcenter+1);
		int rSE=rcenter-1; int zSE=zcenter+1; Cell& southeast = layer.get_cell(rSE,zSE);
		southeast.diagnostic_arrays.area_weighted_charge+=Q*fi*(eta*ksi);
	}
	
	// EAST
	if ( zcenter < layer.z_dim-2 ) {
		fi = 1./(Ncell1*(rcenter+0.5)*layer.n_aver);
		// printf(", to E: %g, area=%g\n", Q*fi*(eta*(1.0-ksi)), (eta*(1.0-ksi)));
		// printf(", rE=%i zE=%i\n", rcenter,zcenter+1);
		int rE=rcenter; int zE=zcenter+1;	Cell& east = layer.get_cell(rE,zE);
		east.diagnostic_arrays.area_weighted_charge+=Q*fi*(eta*(1.0-ksi));
	}
	
	// SOUTH
	if ( rcenter > 0 ) {
		fi = 1./(Ncell1*(rcenter-1+0.5)*layer.n_aver);
		// printf(", to S: %g, area=%g\n", Q*fi*((1.0-eta)*ksi), ((1.0-eta)*ksi));
		// printf(", rS=%i zS=%i\n", rcenter-1,zcenter);
		int rS=rcenter-1; int zS=zcenter;	Cell& south = layer.get_cell(rS,zS);
		south.diagnostic_arrays.area_weighted_charge+=Q*fi*(ksi*(1.0-eta));
	}
#endif
}
