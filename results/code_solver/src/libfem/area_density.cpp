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
#include "grid.h"

#include "fem_solver.h"
#include "debug_printing.h"
#include "fem_debug.h"


	// this is the area weighted density deposition 
	// as an predecessor for the finite element solver method using such area weights
	// and through-cell-wall charge currents
	//
	//									$			
	//						|			$		^	|	
	//						|			$		|	|	
	//	(top_left)|			$		dr|	  (top_right)
	//						|			$		|	|	
	//					A	|B	<-dr->	|	r+1,z+1    C
	// -----------|-----------|-------------- 
	//						|			$	r,z |	r+1,z 				^
	//						|			$			|								|	r_dim
	//~~~~~~~~~~~~|eta~~~~~~~~~~~~~~~~~~~~~~~~|~~~
	//		++++++++++++	$			|	              |
	//		+			D |	 + ksi  E |	           F
	// ---+-------|--+--------|--------------
	//		+		  Q |	 +	$			|	
	//		+			  |	 +	$			|	
	//	 (bottom_left)	$			|	   (bottom_right)
	//		++++++++++++	$			|	
	//		G      	|	H 	$			|	I
	//									$ 
	//							--------> z_dim
	//									$
	//
	// the amount of a charge written to a cell i,j is proportional -
	// here it is, by difference of a sign, exactly the proportion of a standard charge square
	// of dr X dr, where as a cell is dr X dr too
	// additionally, the factor fi for cylinder scaling has to be aaccounted for
	//
	// there are 4 different calculation schemes( top_left, top_right, bottom_left, bottom_right )
	// that have to be applied
	// according to the position of the center of, namely vastly
	// extended, charge surface qu
	//
	// go through all charged particles rather than cells, 
	// so all particles are only touched once - 
	// instead of multiple times for, e.g. neighbouring cells
	//
	// store area weighted charge in cell as diagnostic property


void top_left ( double prtr, double prtz,
		unsigned int rcenter, unsigned int zcenter,
		double Q, unsigned int n_aver ){
#if USE_FEM_SOLVER
	
	// using algorithm delta
	double eta=(zcenter+2)-(prtz-0.5); 
	double ksi=prtr+0.5-rcenter-1;       
	double fi = 1./(Ncell1*(rcenter+0.5)*n_aver);

  GridLayer& layer=global_grid.layers[ELECTRONS];

	iiiprintf ( "using TOP_LEFT;"
						 "assigned to CENTER: %g, area=%g,"
						 " charge %g, fi=%g, eta=%g, ksi=%g,"
						 " to NW: %g, area=%g,"
						 " rNW=%i zNW=%i,"
						 " to W: %g, area=%g,"
						 " rW=%i zW=%i,"
						 " to N: %g, area=%g,"
						 " rN=%i zN=%i\n",
						 Q*fi*((1.0-eta)*(1.0-ksi)), (1.0-eta)*(1.0-ksi),
						 Q, fi, eta, ksi,
						 Q*fi*(eta*ksi), (eta*ksi),
						 rcenter+1, zcenter-1,
						 Q*fi*(eta*(1.0-ksi)), (eta*(1.0-ksi)),
						 rcenter, zcenter-1,
						 Q*fi*(ksi*(1.0-eta)), ((1.0-eta)*ksi),
						 rcenter+1, zcenter);

	// get relevant neighbouring cells 
	// CENTER, NORTH, WEST, NORTHWEST
	// calculating area charge for four cells successively
	
	// CENTER
	Cell& middle = layer.get_cell(rcenter,zcenter);
	middle.fem.area_weighted_charge+=Q*fi*((1.0-ksi)*(1.0-eta));

	// NORTHWEST
	if ( (rcenter < global_grid.r_dim-2) && (zcenter > 0) ) {
		fi = 1./(Ncell1*(rcenter+1+0.5)*n_aver);
		int rNW=rcenter+1; int zNW=zcenter-1; Cell& northwest = layer.get_cell(rNW,zNW);
		northwest.fem.area_weighted_charge+=Q*fi*(eta*ksi);
	}
	// WEST
	if ( zcenter > 0 ) {
		fi = 1./(Ncell1*(rcenter+0.5)*n_aver);
		int rW=rcenter; int zW=zcenter-1;	Cell& west = layer.get_cell(rW,zW);
		west.fem.area_weighted_charge+=Q*fi*(eta*(1.0-ksi));
	}
	// NORTH
	if ( rcenter < global_grid.r_dim-2 ) {
		fi = 1./(Ncell1*(rcenter+1+0.5)*n_aver);
		int rN=rcenter+1; int zN=zcenter;	Cell& north = layer.get_cell(rN,zN);
		north.fem.area_weighted_charge+=Q*fi*(ksi*(1.0-eta));
	}

#endif
}

void top_right (double prtr, double prtz,
		unsigned int rcenter, unsigned int zcenter,
		double Q, unsigned int n_aver  ){
#if USE_FEM_SOLVER
	
	// using algorithm gamma
	double eta=prtz+0.5-zcenter-1;
	double ksi=prtr+0.5-rcenter-1;
	double fi = 1./(Ncell1*(rcenter+0.5)*n_aver);

  GridLayer& layer=global_grid.layers[ELECTRONS];

	iiprintf ( "using TOP_RIGHT; "
	           "assigned to CENTER: %g, area=%g,"
	           " charge %g, fi=%g, eta=%g, ksi=%g,"
	           " to NE: %g, area=%g,"
	           " rNE=%i zNE=%i,"
	           " to E: %g, area=%g,"
	           " rE=%i zE=%i,"
	           " to N: %g, area=%g,"
	           " rN=%i zN=%i\n",
             Q*fi*((1.0-eta)*(1.0-ksi)), (1.0-eta)*(1.0-ksi),
             Q, fi, eta, ksi,
             Q*fi*(eta*ksi), (eta*ksi),
             rcenter+1, zcenter+1,
             Q*fi*(eta*(1.0-ksi)), (eta*(1.0-ksi)),
             rcenter, zcenter+1,
             Q*fi*((1.0-eta)*ksi), ((1.0-eta)*ksi),
             rcenter+1, zcenter);

	// get relevant neighbouring cells 
	// NORTH, WEST, NORTHWEST
	// calculating area charge for four cells successively
	
	// CENTER
	Cell& middle = layer.get_cell(rcenter,zcenter);
	middle.fem.area_weighted_charge+=Q*fi*((1.0-ksi)*(1.0-eta));
	
	// NORTHEAST
	if ( (rcenter < global_grid.r_dim-2) && (zcenter < global_grid.z_dim-2) ) {
		fi = 1./(Ncell1*(rcenter+1+0.5)*n_aver);
		int rNE=rcenter+1; int zNE=zcenter+1; Cell& northeast = layer.get_cell(rNE,zNE);
		northeast.fem.area_weighted_charge+=Q*fi*(eta*ksi);
	}
	// EAST
	if ( zcenter < global_grid.z_dim-2 ) {
		fi = 1./(Ncell1*(rcenter+0.5)*n_aver);
		int rE=rcenter; int zE=zcenter+1;	Cell& east = layer.get_cell(rE,zE);
		east.fem.area_weighted_charge+=Q*fi*(eta*(1.0-ksi));
	}
	// NORTH
	if ( rcenter < global_grid.r_dim-2 ) {
		fi = 1./(Ncell1*(rcenter+0.5)*n_aver);
		int rN=rcenter+1; int zN=zcenter;	Cell& north = layer.get_cell(rN,zN);
		north.fem.area_weighted_charge+=Q*fi*(ksi*(1.0-eta));
	}

#endif
}

void bottom_left ( double prtr, double prtz,
		unsigned int rcenter, unsigned int zcenter,
		double Q, unsigned int n_aver  ){
#if USE_FEM_SOLVER
	
	// using algorithm beta
	double ksi=(rcenter-1+1)-(prtr-0.5);
	double eta=zcenter-1+1-(prtz-0.5);
	double fi = 1./(Ncell1*(rcenter+0.5)*n_aver);

  GridLayer& layer=global_grid.layers[ELECTRONS];

	iiiprintf ( "using BOTTOM_LEFT; "
	            "assigned to CENTER: %g, area=%g"
	            " charge %g, fi=%g, eta=%g, ksi=%g"
	            " to SW: %g, area=%g"
	            " rS=%i zS=%i"
	            " to S: %g, area=%g"
	            " rS=%i zS=%i"
	            " to W: %g, area=%g"
	            " rW=%i zW=%i\n",
					    Q*fi*((1.0-eta)*(1.0-ksi)), (1.0-eta)*(1.0-ksi),
					    Q, fi, eta, ksi,
					    Q*fi*(eta*ksi), (eta*ksi),
					    rcenter-1, zcenter-1,
					    Q*fi*((1.0-eta)*ksi), ((1.0-eta)*ksi),
					    rcenter-1,zcenter,
					    Q*fi*(eta*(1.0-ksi)), (eta*(1.0-ksi)),
					    rcenter, zcenter-1);

	// get relevant neighbouring cells 
	// WEST, SOUTH, SOUTH-WEST
	// calculating area charge for four cells successively
	
	// CENTER
	Cell& middle = layer.get_cell(rcenter,zcenter);
	middle.fem.area_weighted_charge+=Q*fi*((1.0-ksi)*(1.0-eta));

	// SOUTH-WEST
	if ( (rcenter > 0) && (zcenter > 0) ) {
		fi = 1./(Ncell1*(rcenter-1+0.5)*n_aver);
		int rSW=rcenter-1; int zSW=zcenter-1; Cell& southwest = layer.get_cell(rSW,zSW);
		southwest.fem.area_weighted_charge+=Q*fi*(eta*ksi);
	}
	// WEST
	if ( zcenter > 0 ) {
		fi = 1./(Ncell1*(rcenter+0.5)*n_aver);
		int rW=rcenter; int zW=zcenter-1; Cell& west = layer.get_cell(rW,zW);
		west.fem.area_weighted_charge+=Q*fi*(eta*(1.0-ksi));
	}
	// SOUTH
	if ( rcenter > 0 ) {
		fi = 1./(Ncell1*(rcenter-1+0.5)*n_aver);
		int rS=rcenter-1; int zS=zcenter; Cell& south = layer.get_cell(rS,zS);
		south.fem.area_weighted_charge+=Q*fi*(ksi*(1.0-eta));
	}

#endif
}

void bottom_right ( double prtr, double prtz,
		unsigned int rcenter, unsigned int zcenter,
		double Q, unsigned int n_aver  ){
#if USE_FEM_SOLVER
	
	// using algorithm alpha
	double eta=prtz+0.5-zcenter-1;
	double ksi=(rcenter-1+1)-(prtr-0.5);
	double fi = 1./(Ncell1*(rcenter+0.5)*n_aver);

  GridLayer& layer=global_grid.layers[ELECTRONS];

  iiiprintf ( "using BOTTOM_RIGHT; "
	            "assigned to CENTER: %g, area=%g"
	            " charge %g, fi=%g, eta=%g, ksi=%g"
	            " to SE: %g, area=%g"
	            " rSE=%i zSE=%i"
	            " to E: %g, area=%g"
	            " rE=%i zE=%i"
	            " to S: %g, area=%g"
	            " rS=%i zS=%i\n",
							Q*fi*((1.0-eta)*(1.0-ksi)), (1.0-eta)*(1.0-ksi),
							Q, fi, eta, ksi,
							Q*fi*(eta*ksi), (eta*ksi),
							rcenter-1, zcenter+1,
							Q*fi*(eta*(1.0-ksi)), (eta*(1.0-ksi)),
							rcenter, zcenter+1,
							Q*fi*((1.0-eta)*ksi), ((1.0-eta)*ksi),
							rcenter-1, zcenter);

 	// get relevant neighbouring cells 
	// calculating area charge for four cells successively
	// EAST, SOUTH, SOUTH-EAST
	
	// CENTER
	Cell& middle = layer.get_cell(rcenter,zcenter);
	middle.fem.area_weighted_charge+=Q*fi*((1.0-ksi)*(1.0-eta));
	// SOUTH-EAST
	if ( (rcenter > 0) && (zcenter < global_grid.z_dim-2) ) {
		fi = 1./(Ncell1*(rcenter-1+0.5)*n_aver);
		int rSE=rcenter-1; int zSE=zcenter+1; Cell& southeast = layer.get_cell(rSE,zSE);
		southeast.fem.area_weighted_charge+=Q*fi*(eta*ksi);
	}
	// EAST
	if ( zcenter < global_grid.z_dim-2 ) {
		fi = 1./(Ncell1*(rcenter+0.5)*n_aver);
		int rE=rcenter; int zE=zcenter+1;	Cell& east = layer.get_cell(rE,zE);
		east.fem.area_weighted_charge+=Q*fi*(eta*(1.0-ksi));
	}
	// SOUTH
	if ( rcenter > 0 ) {
		fi = 1./(Ncell1*(rcenter-1+0.5)*n_aver);
		int rS=rcenter-1; int zS=zcenter;	Cell& south = layer.get_cell(rS,zS);
		south.fem.area_weighted_charge+=Q*fi*(ksi*(1.0-eta));
	}

#endif
}

void area_density(){
#if USE_FEM_SOLVER

	double eta,ksi,fi,qu;
	unsigned int rO,zO,n_aver;

	for (auto& layer: global_grid.layers){			//LAYERS, EXCEPT NEUTRALS
		if ( layer.name=="neutrals" ){ } else {		//ONLY CHARGED PARTICLES
			
			//CHARGE IN UNITS OF E0
			if ( layer.name=="electrons"		 ) { qu=-1.0;
			} else if ( layer.name=="ions"	 ) { qu= 1.0; 
			} else if ( layer.name=="nions"  ) { qu=-1.0;
			} else if ( layer.name=="ions2p" ) { qu= 2.0;
			}

		  if ( layer.n_aver < 1) { layer.n_aver = 1; }
        n_aver = layer.n_aver;
      
			  for ( unsigned int rG = 0; rG < layer.r_dim; ++rG ){

					for ( unsigned int zG = 0; zG < layer.z_dim; ++zG ){

						// mother cell
						Cell& center = layer.get_cell(rG,zG);

						for ( unsigned int ptid = 0; ptid < center.size(); ++ptid ){

							// get particles from cell
							Particle& pt=center.particles[ptid];

							double ptr=pt.r; double ptz=pt.z;
							double rest=ptr-rG, zest=ptz-zG;
				
							iiiprintf ( "\n>> cell & algorithm info: "
													"r=%i/(%i), z=%i/(%i); center pos. R=%g, Z=%g; "
													"ptid=%i/(%i), "
													"particle pos. r=%g, z=%g,"
													" cell pos. rest=%g, zest=%g\n",
													rG, layer.r_dim, layer.z_dim, zG, rG+0.5, zG+0.5, ptid, center.size(),
													ptr, ptz, rest, zest);


							// ALGORITHMS 
							if ( ptr < rG+0.5 ){ 
								if ( ptz < zG+0.5 ){	 // BOTTOM_LEFT
									bottom_left (ptr,ptz,rG,zG,qu,n_aver);
	    		    	} else {							 // BOTTOM_RIGHT
									bottom_right (ptr,ptz,rG,zG,qu,n_aver);
	    		      }
	    		    } else {								
					  		if ( pt.z < zG+0.5 ){  // TOP_LEFT
									top_left (ptr,ptz,rG,zG,qu,n_aver);
	    		      } else {							 // TOP_RIGHT
									top_right (ptr,ptz,rG,zG,qu,n_aver);
					  	  }
					  	} 

						} // particle ids
	 				} // zdim
				} // rdim
			} // except neutrals
		} // layers

#endif
}
