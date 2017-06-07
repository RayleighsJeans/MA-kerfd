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
  std::vector<std::vector<double>> areaw_charge = { {0.0} };
  areaw_charge.resize(layer.r_dim);
  for ( unsigned int r = 0; r < layer.r_dim; ++r){
    areaw_charge[r].resize(layer.z_dim);
    for (unsigned int z = 0; z < layer.z_dim; ++z){
      Cell& cell = layer.get_cell(r,z);
      areaw_charge[r][z] = cell.fem.area_weighted_charge;
    }
  }
  reduce(areaw_charge);

	int mpi_rank=get_rank();
  if (mpi_rank==0){
    FILE  *file = fopen(dat_name.c_str(), "w");
    if ( !file ){ 
      printf("could not open %s\n",dat_name.c_str()); 
      exit(EXIT_FAILURE);
    }

    for ( unsigned int r = 0; r < layer.r_dim; ++r){
      for (unsigned int z = 0; z < layer.z_dim; ++z){
        fprintf(file, "% .5e ", areaw_charge[r][z]);
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
  std::vector<std::vector<double>> current = {{ 0.0 }};
  current.resize(layer.r_dim);
  for ( unsigned int r = 0; r < layer.r_dim; ++r){
    current[r].resize(layer.z_dim);
    for (unsigned int z = 0; z < layer.z_dim; ++z){
		  Cell& cell = layer.get_cell(r,z);
      current[r][z] += cell.fem.current_bottom +
                       cell.fem.current_left   -
                       cell.fem.current_top    -
                       cell.fem.current_right;
      if ( z > 0 ) {
        iiiprintf ( ">> out cell current: leftcell, current=%g, r=%d, z-1=%d\n",
                    current[r][z], r, z-1 );
        Cell& leftcell = layer.get_cell(r,z-1);
        current[r][z] += leftcell.fem.current_right;
      } else if ( r > 0 ) {
        iiiprintf ( ">> out cell current: bottomcell, current=%g, r-1=%d, z=%d\n",
                    current[r][z], r-1, z );
        Cell& bottomcell = layer.get_cell(r-1,z);
        current[r][z] += bottomcell.fem.current_top;
      } else if ( r < layer.r_dim-1 ) {
        iiiprintf ( ">> out cell current: topcell, current=%g, r+1=%d, z=%d\n",
                    current[r][z], r+1, z );
        Cell& topcell = layer.get_cell(r+1,z);
        current[r][z] -= topcell.fem.current_bottom;
      } else if ( z < layer.z_dim-1 ) {
        iiiprintf ( ">> out cell current: rightcell, current=%g, r=%d, z+1=%d\n",
                   current[r][z], r, z+1 );
        Cell& rightcell = layer.get_cell(r,z+1);
        current[r][z] -= rightcell.fem.current_left;
      }
      iiiprintf ( ">> out cell current: current=%g, r=%d, z=%d\n",
                 current[r][z], r, z );
    }
  }
  reduce(current);

  int mpi_rank=get_rank();
  if (mpi_rank==0){
    FILE  *file = fopen(dat_name.c_str(), "w");
    if ( !file ){ 
      printf("could not open %s\n",dat_name.c_str()); 
      exit(EXIT_FAILURE);
    }
    for ( unsigned int r = 0; r < layer.r_dim; ++r){
      for (unsigned int z = 0; z < layer.z_dim; ++z){
        fprintf(file, " % .5e",	current[r][z]);
      }
	  	fprintf(file, "\n");
    }
   fclose(file);
  }

#endif
}

void out_efield ( std::string dat_name1, std::string dat_name2 ) {
#if USE_FEM_SOLVER

/**********************************************************************************/
/*********************** MY NEW ELECTRIC FIELD ************************************/
/**********************************************************************************/

  GridLayer& layer=global_grid.layers[ELECTRONS];
  std::vector<std::vector<double>> efield = {{0.0}}; 
  efield.resize(layer.r_dim);
  for ( unsigned int r = 0; r < layer.r_dim; ++r){ 
    efield[r].resize(layer.z_dim);
    for (unsigned int z = 0; z < layer.z_dim; ++z){
      Cell& cell = layer.get_cell(r,z);
      efield[r][z] += cell.fem.efield_left   +
                      cell.fem.efield_bottom -
                      cell.fem.efield_top    -
                      cell.fem.efield_right;
      iiiprintf(">> out efield: centercell, efield=%g, r=%d/%d, z=%d/%d\n",
                efield[r][z], r, layer.r_dim-1, z, layer.z_dim-1);
      if ( z > 0 ) {
       iiiprintf(">> out efield: leftcell, efield=%g, r=%d/%d, z=%d/%d\n",
                  efield[r][z], r, layer.r_dim-1, z-1, layer.z_dim-1);
        Cell& leftcell = layer.get_cell(r,z-1);
        efield[r][z] += leftcell.fem.efield_right;
      } 
      if ( r > 0 ) {
        iiiprintf(">> out efield: bottomcell, efield=%g, r=%d/%d, z=%d/%d\n",
                  efield[r][z], r-1, layer.r_dim-1, z, layer.z_dim-1);
        Cell& bottomcell = layer.get_cell(r-1,z);
        efield[r][z] += bottomcell.fem.efield_top;
      } 
      if ( r < layer.r_dim-1 ) {
        iiiprintf(">> out efield: topcell, efield=%g, r=%d/%d, z=%d/%d\n",
                  efield[r][z], r+1, layer.r_dim-1, z, layer.z_dim-1);
        Cell& topcell = layer.get_cell(r+1,z);
        efield[r][z] -= topcell.fem.efield_bottom;
      }
      if ( z < layer.z_dim-1 ) {
        iiiprintf(">> out efield: rightcell, efield=%g, r=%d/%d, z=%d/%d\n",
                  efield[r][z], r, layer.r_dim-1, z+1, layer.z_dim-1);
        Cell& rightcell = layer.get_cell(r,z+1);
        efield[r][z] -= rightcell.fem.efield_left;
      }
    }
  }
  reduce(efield);
  
  int mpi_rank=get_rank();
  if (mpi_rank==0){
    FILE  *file1 = fopen(dat_name1.c_str(), "w");
    if ( !file1 ){ 
      printf("could not open %s\n",dat_name1.c_str()); 
      exit(EXIT_FAILURE); }
      
    for ( unsigned int r = 0; r < layer.r_dim; ++r){ 
        for (unsigned int z = 0; z < layer.z_dim; ++z){
          fprintf(file1, " % .5e",	efield[r][z]);
        }
        fprintf(file1, "\n");
      }
    fclose(file1);

/**********************************************************************************/
/**************** NORMAL ELECTRIC FIELD OUTPUT FOR COMPARISON *********************/
/**********************************************************************************/

    FILE  *file2 = fopen(dat_name2.c_str(), "w");
    if ( !file2 ){ 
      printf("could not open %s\n",dat_name2.c_str()); 
      exit(EXIT_FAILURE);
    }

    int i, j;
    int NZ = global_grid.mesh_z_dim;
    double field=0.0;
    for (i = 1; i < layer.r_dim; ++i) {
      for (j = 1; j < layer.z_dim; ++j) {
  
        field += E_grid[i*NZ + j].r; // left
        field += E_grid[i*NZ + j].z; // bottom
        
        if ( j > 1 ) {
          field += E_grid[i*NZ + (j-1)].z; // left grid point
        } else if ( i > 1 ) {
          field += E_grid[(i-1)*NZ + j].r; // bottom grid point
        } else if ( i < layer.r_dim-1 ) {
          field -= E_grid[(i+1)*NZ + j].r; // top grid point
        } else if ( j < layer.z_dim-1 ) {
          field -= E_grid[i*NZ + (j+1)].z; // right grid point
        }
  
        fprintf(file2, " % .5e",	field);
        field=0; // RESET VALUE
      }
      fprintf(file2, "\n");
    }
    fclose(file2);
    field=0.0; // CLEAR
  }

#endif
}

