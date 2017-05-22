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
#include "grid.h"
#include "mpi_wrapper.h"

#include "fem_solver.h"
#include "debug_printing.h"
#include "fem_debug.h"

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
        fprintf(file, "% .5e ", cell.fem.area_weighted_charge);               
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
        current += cell.fem.current_bottom +
                   cell.fem.current_left   -
                   cell.fem.current_top    -
                   cell.fem.current_right;
        
        if ( z > 0 ) {
          iiprintf ( ">> out cell current: leftcell, current=%g, r=%d, z-1=%d\n",
                     current, r, z-1 );
          Cell& leftcell = layer.get_cell(r,z-1);
          current += leftcell.fem.current_right;
        } else if ( r > 0 ) {
          iiprintf ( ">> out cell current: bottomcell, current=%g, r-1=%d, z=%d\n",
                     current, r-1, z );
          Cell& bottomcell = layer.get_cell(r-1,z);
          current += bottomcell.fem.current_top;
        } else if ( r < layer.r_dim-1 ) {
          iiprintf ( ">> out cell current: topcell, current=%g, r+1=%d, z=%d\n",
                     current, r+1, z );
          Cell& topcell = layer.get_cell(r+1,z);
          current -= topcell.fem.current_bottom;
        } else if ( z < layer.z_dim-1 ) {
          iiprintf ( ">> out cell current: rightcell, current=%g, r=%d, z+1=%d\n",
                     current, r, z+1 );
          Cell& rightcell = layer.get_cell(r,z+1);
          current -= rightcell.fem.current_left;
        }
      
      iiprintf ( ">> out cell current: current=%g, r=%d, z=%d\n",
                 current, r, z );
      fprintf(file, " % .5e",	current);
      current=0; // RESET VALUE
      }
			fprintf(file, "\n");
    }
    fclose(file);
  }

#endif
}

void out_efield ( std::string dat_name1, std::string dat_name2 ) {
#if USE_FEM_SOLVER

  GridLayer& layer=global_grid.layers[ELECTRONS];
  int mpi_rank=get_rank();
  if (mpi_rank==0){
    FILE  *file1 = fopen(dat_name1.c_str(), "w");
    if ( !file1 ){ 
      printf("could not open %s\n",dat_name1.c_str()); 
      exit(EXIT_FAILURE);
    }

    double efield;
    for ( unsigned int r = 0; r < layer.r_dim; ++r){
      for (unsigned int z = 0; z < layer.z_dim; ++z){
        Cell& cell = layer.get_cell(r,z);
        efield  += cell.fem.efield_left   +
                   cell.fem.efield_bottom -
                   cell.fem.efield_top    -
                   cell.fem.efield_right;
        
        if ( z > 0 ) {
          iiprintf ( ">> out efield: leftcell, efield=%g, r=%d, z-1=%d\n",
                     efield, r, z-1 );
          Cell& leftcell = layer.get_cell(r,z-1);
          efield += leftcell.fem.efield_right;
        } else if ( r > 0 ) {
          iiprintf ( ">> out efield: bottomcell, efield=%g, r-1=%d, z=%d\n",
                     efield, r-1, z );
          Cell& bottomcell = layer.get_cell(r-1,z);
          efield += bottomcell.fem.efield_top;
        } else if ( r < layer.r_dim-1 ) {
          iiprintf ( ">> out efield: topcell, efield=%g, r+1=%d, z=%d\n",
                     efield, r+1, z );
          Cell& topcell = layer.get_cell(r+1,z);
          efield -= topcell.fem.efield_bottom;
        } else if ( z < layer.z_dim-1 ) {
          iiprintf ( ">> out efield: rightcell, efield=%g, r=%d, z+1=%d\n",
                     efield, r, z+1 );
          Cell& rightcell = layer.get_cell(r,z+1);
          efield -= rightcell.fem.efield_left;
        }
      
      iiprintf ( ">> out efield: efield=%g, r=%d, z=%d\n",
                 efield, r, z );

      fprintf(file1, " % .5e",	efield);;
      efield=0; // RESET VALUE
      }
      fprintf(file1, "\n");
    }
    fclose(file1);
    efield=0; // CLEAR

   FILE  *file2 = fopen(dat_name2.c_str(), "w");
    if ( !file2 ){ 
      printf("could not open %s\n",dat_name2.c_str()); 
      exit(EXIT_FAILURE);
    }

    int i, j;
    int NZ = global_grid.mesh_z_dim;
  
    for (i = 1; i < layer.r_dim; ++i) {
      for (j = 1; j < layer.z_dim; ++j) {
  
        efield += E_grid[i*NZ + j].r; // left
        efield += E_grid[i*NZ + j].z; // bottom
        
        if ( j > 1 ) {
          efield += E_grid[i*NZ + (j-1)].z; // left grid point
        } else if ( i > 1 ) {
          efield += E_grid[(i-1)*NZ + j].r; // bottom grid point
        } else if ( i < layer.r_dim-1 ) {
          efield -= E_grid[(i+1)*NZ + j].r; // top grid point
        } else if ( j < layer.z_dim-1 ) {
          efield -= E_grid[i*NZ + (j+1)].z; // right grid point
        }
  
        fprintf(file2, " % .5e",	efield);;
        efield=0; // RESET VALUE
      }
      fprintf(file2, "\n");
    }
    fclose(file2);
    efield=0; // CLEAR
  }

#endif
}

