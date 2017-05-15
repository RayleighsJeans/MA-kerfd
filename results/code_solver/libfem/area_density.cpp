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

void area_density(){
#if USE_FEM_SOLVER

	// this is THE TEST for the area weighted density deposition 
	// as an predecessor for the finite element solver method using such area weights
	// and through-cell-wall charge currents
	//
	//								$			
	//					|			$		^	|	
	//					|			$		|	|	
	//(top_left)|			$		dr|	  (top_right)
	//					|			$		|	|	
	//				A	|B	<-dr->	|	r+1,z+1    C
	// ---------|-----------|-------------- 
	//					|			$	r,z |	r+1,z 				^
	//					|			$			|								|	r_dim
	//~~~~~~~~~~|eta~~~~~~~~~~~~~~~~~~~~~~~~|~~~
	//	++++++++++++	$			|	              |
	//	+			D |	 + ksi  E |	           F
	// -+-------|--+--------|--------------
	//	+		  Q |	 +	$			|	
	//	+			  |	 +	$			|	
	// (bottom_left)	$			|	   (bottom_right)
	//	++++++++++++	$			|	
	//	G      	|	H 	$			|	I
	//								$ 
	//						--------> z_dim
	//								$
	
	double eta,ksi,fi,qu;
	unsigned int rO,zO;

	// the amount of a charge written to a cell i,j is proportional -
	// here it is, by difference of a sign, exactly the proportion of a standard charge square
	// of dr X dr, where as a cell is dr X dr too
	// additionally, the factor fi for cylinder scaling has to be aaccounted for
	//
	// there are 4 different calculation schemes( top_left, top_right, bottom_left, bottom_right )
	// that have to be applied
	// according to the position of the center of, namely vastly
	// extended, charge surface Q
	//
	// go through all charged particles rather than cells, 
	// so all particles are only touched once - 
	// instead of multiple times for, e.g. neighbouring cells
	//
	// store area weighted charge in cell as diagnostic property

	for (auto& layer: global_grid.layers){ //LAYERS, EXCEPT NEUTRALS
		if ( layer.name=="neutrals" ){ } else { //ONLY CHARGED PARTICLES
			//CHARGE IN UNITS OF E0
			if ( layer.name=="electrons" ){ qu=-1.0; } else if ( layer.name=="ions" ) { qu=1.0; 
			} else if ( layer.name=="nions" ) { qu=-1.0; } else if ( layer.name=="ions2p" ) { qu=2.0; } 
		  if ( layer.n_aver < 1) { layer.n_aver = 1; }
			  for ( unsigned int rG = 0; rG < layer.r_dim; ++rG ){
					for ( unsigned int zG = 0; zG < layer.z_dim; ++zG ){
						
						// printf(">> cell & algorithm info: ");
						// mother cell
						Cell& center = layer.get_cell(rG,zG);
						// printf("r=%i,z=%i; center pos. R=%g, Z=%g; \n", rG, zG, rG+0.5, zG+0.5);
						for ( unsigned int ptid = 0; ptid < center.size(); ++ptid ){
							// get particles from cell
							Particle& pt=center.particles[ptid];
							double ptr=pt.r; double ptz=pt.z;
							// printf("ptid=%i/(%i), \n", ptid, center.size());
							// printf("particle pos. r=%g, z=%g; \n", pt.r, pt.z);

							// ALGORITHMS 
							if ( ptr < rG+0.5 ){ 
								if ( ptz < zG+0.5 ){	 // BOTTOM_LEFT
									bottom_left (layer,ptr,ptz,rG,zG,qu);
	    		    	} else {							 // BOTTOM_RIGHT
									bottom_right (layer,ptr,ptz,rG,zG,qu);
	    		      }
	    		    } else {								
					  		if ( pt.z < zG+0.5 ){  // TOP_LEFT
									top_left (layer,ptr,ptz,rG,zG,qu);
	    		      } else {							 // TOP_RIGHT
									top_right (layer,ptr,ptz,rG,zG,qu);
					  	  }
					  	} // printf(" \n");
						} // particle ids
	 				} // zdim
				} // rdim
			} // except neutrals
		} // layers
#endif
}
