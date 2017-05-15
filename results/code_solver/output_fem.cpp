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

void out_cellcurrent ( std::string dat_name ) {
#if USE_FEM_SOLVER

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
				current = ( cell.diagnostic_arrays.current_top +
										cell.diagnostic_arrays.current_bottom +
									  cell.diagnostic_arrays.current_right +
										cell.diagnostic_arrays.current_left )/4;
				fprintf(file, " % .5e",	current);
      }
			fprintf(file, "\n");
    }
    fclose(file);
  }

#endif
}
