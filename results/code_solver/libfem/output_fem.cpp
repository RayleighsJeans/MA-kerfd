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

void out_area_density(std::string dat_name){
#if USE_FEM_SOLVER

	GridLayer& layer=global_grid.layers[ELECTRONS];
	int mpi_rank=get_rank();
  if (mpi_rank==0){
    FILE  *file = fopen(dat_name.c_str(), "w");
    if ( !file ){ 
      printf("could not open %s\n",dat_name.c_str()); 
      exit(EXIT_FAILURE);
    }

    for ( unsigned int r = 0; r < layer.r_dim; ++r){
      for (unsigned int z = 0; z < layer.z_dim; ++z){
			  Cell& cell = layer.get_cell(r,z);
        fprintf(file, "% .5e ", cell.diagnostic_arrays.area_weighted_charge);               
      }
      fprintf(file, "\n");        
    }
    fclose(file);
  }
#endif
}

void out_cellcurrent(std::string dat_name){
	GridLayer& layer=global_grid.layers[ELECTRONS];
	int mpi_rank=get_rank();
  if (mpi_rank==0){
    FILE  *file = fopen(dat_name.c_str(), "w");
    if ( !file ){ 
      printf("could not open %s\n",dat_name.c_str()); 
      exit(EXIT_FAILURE);
    }

		double current;
    for ( unsigned int r = 0; r < layer.r_dim; ++r){
      for (unsigned int z = 0; z < layer.z_dim; ++z){
			  Cell& cell = layer.get_cell(r,z);
				current = ( cell.diagnostic_arrays.current_north +
										cell.diagnostic_arrays.current_south +
									  cell.diagnostic_arrays.current_east +
										cell.diagnostic_arrays.current_west )/4;

				fprintf(file, " % .5e",	current);

		//		if ( r == 0 ) {
		//			fprintf(file, " % 0.00000e0 % .5e",
		//					cell.diagnostic_arrays.current_south);
		//		} else {
		//			fprintf(file, "% .5e  % .5e",
		//					cell.diagnostic_arrays.current_west,
		//					cell.diagnostic_arrays.current_east);
		//		}

      }
      
		//	if (r!=0) { fprintf(file, "\n"); }
    //  for (unsigned int z2 = 0; z2 < layer.z_dim; ++z2) {
		//		Cell& cell2 = layer.get_cell(r,z2);
		//		if ( (r != layer.r_dim-1) && (r != 0) ) {
		//			Cell& topcell = layer.get_cell(r+1,z2);
		//			fprintf(file, "% .5e  % .5e",
		//					topcell.diagnostic_arrays.current_south,
		//					cell2.diagnostic_arrays.current_north);
		//		} else if (r!=0) {
		//			fprintf(file, "% .5e % 0.00000e0",
		//					cell2.diagnostic_arrays.current_north); }
		//	}

			fprintf(file, "\n");
    }
    fclose(file);
  }
}
