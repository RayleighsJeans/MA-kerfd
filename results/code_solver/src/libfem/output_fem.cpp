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
//#define DEBUG_LEVEL DEBUG_ERROR

void output_fem ( Field E_grid[] ) {
  const int MAX_CHAR_SIZE = 256;
  char numb[MAX_CHAR_SIZE];

  snprintf(numb, MAX_CHAR_SIZE, "%08i", nstep);
  std::string path = "out/";

  auto get_full_path = [&numb,&path]( std::string filename ){ 
    std::string storage;
    storage = path + filename + std::string(numb) + ".dat";
    return storage; 
  };

  out_area_density  ( get_full_path("aw_charge_") );
  out_cellcurrent   ( get_full_path("cell_current_") );
  out_efield        ( get_full_path("fem_efield_"),
                      get_full_path("stand_efield_"),
                      get_full_path("diff_fields_"),
                      E_grid,
                      get_full_path("old_field_") );
  out_vector_femfield ( get_full_path("vector_femfield_") );
  out_vector_stdfield ( E_grid, get_full_path("vector_stdfield_") );

}

void out_area_density( std::string dat_name ){
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
  if (mpi_rank==0)
  { FILE  *file = fopen(dat_name.c_str(), "w");
    if ( !file )
    { printf("could not open %s\n",dat_name.c_str()); 
      exit(EXIT_FAILURE); }

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
      current[r][z] = cell.fem.current_bottom +
                      cell.fem.current_left   -
                      cell.fem.current_top    -
                      cell.fem.current_right;
      if ( z > 0 ) {
        iiiprintf ( ">> out cell current: leftcell, "
                    "current=%g, r=%d, z-1=%d\n",
                    current[r][z], r, z-1 );
        Cell& leftcell = layer.get_cell(r,z-1);
        current[r][z] += leftcell.fem.current_right;
      } else if ( r > 0 ) {
        iiiprintf ( ">> out cell current: bottomcell, "
                    "current=%g, r-1=%d, z=%d\n",
                    current[r][z], r-1, z );
        Cell& bottomcell = layer.get_cell(r-1,z);
        current[r][z] += bottomcell.fem.current_top;
      } else if ( r < layer.r_dim-1 ) {
        iiiprintf ( ">> out cell current: topcell, "
                    "current=%g, r+1=%d, z=%d\n",
                    current[r][z], r+1, z );
        Cell& topcell = layer.get_cell(r+1,z);
        current[r][z] -= topcell.fem.current_bottom;
      } else if ( z < layer.z_dim-1 ) {
        iiiprintf ( ">> out cell current: rightcell, "
                    "current=%g, r=%d, z+1=%d\n",
                    current[r][z], r, z+1 );
        Cell& rightcell = layer.get_cell(r,z+1);
        current[r][z] -= rightcell.fem.current_left;
      }
      iiiprintf ( ">> out cell current: current=%g,"
                  " r=%d, z=%d\n", current[r][z], r, z );
    }
  }
  reduce(current);

  int mpi_rank=get_rank();
  if (mpi_rank==0)
  { FILE  *file = fopen(dat_name.c_str(), "w");
    if ( !file )
    {  printf("could not open %s\n",dat_name.c_str()); 
      exit(EXIT_FAILURE); }
    
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

void out_efield ( std::string dat_name1,
                  std::string dat_name2,
                  std::string dat_name3,
                  Field E_grid[],
                  std::string dat_name4 ) {
#if USE_FEM_SOLVER
  GridLayer& layer=global_grid.layers[ELECTRONS];
  std::vector<std::vector<double>> efield = {{0.0}}; 
  std::vector<std::vector<double>> oldfield = {{0.0}}; 
  efield.resize(layer.r_dim);
  oldfield.resize(layer.r_dim);
  for ( unsigned int r = 0; r < layer.r_dim; ++r){ 
    efield[r].resize(layer.z_dim);
    oldfield[r].resize(layer.z_dim);
    for (unsigned int z = 0; z < layer.z_dim; ++z){
      Cell& cell = layer.get_cell(r,z);
      efield[r][z] = ( cell.fem.efield_top + cell.fem.efield_left -
                       cell.fem.efield_bottom - cell.fem.efield_right )/dr;
      oldfield[r][z] = ( cell.fem.oldfield_top + cell.fem.oldfield_left -
                         cell.fem.oldfield_bottom - cell.fem.oldfield_right )/dr;
    }
  }
  reduce(efield);
  reduce(oldfield);
  
  /************ ELECTRIC FIELD OUTPUT ****************/
  int mpi_rank=get_rank();
  if (mpi_rank==0)
  { FILE  *file1 = fopen(dat_name1.c_str(), "w");
    if ( !file1 ) 
    { printf("could not open %s\n",dat_name1.c_str()); 
      exit(EXIT_FAILURE); }
      
    for ( unsigned int r = 0; r < layer.r_dim; ++r){ 
        for (unsigned int z = 0; z < layer.z_dim; ++z){
          fprintf(file1, " % .5e",	efield[r][z]);
        }
        fprintf(file1, "\n");
      }
    fclose(file1);
    
    /********** OLD ELECTRIC FIELD OUTPUT ************/
    FILE  *file4 = fopen(dat_name4.c_str(), "w");
      if ( !file4 ) 
      { printf("could not open %s\n",dat_name4.c_str()); 
        exit(EXIT_FAILURE); }
        
      for ( unsigned int r = 0; r < layer.r_dim; ++r){ 
          for (unsigned int z = 0; z < layer.z_dim; ++z){
            fprintf(file4, " % .5e",	oldfield[r][z]);
          }
          fprintf(file4, "\n");
        }
      fclose(file4);

    std::vector<std::vector<double>> diff_fields = {{ 0.0 }};
    std::vector<std::vector<double>> std_efield = {{ 0.0 }};
  
    diff_fields.resize(layer.r_dim);
    std_efield.resize(layer.r_dim);
    for ( unsigned int r = 0; r < layer.r_dim; ++r)
    { std_efield[r].resize(layer.z_dim);
      diff_fields[r].resize(layer.z_dim);}

    /************ get old field from funct ***********/
    std_efield = out_std_efield ( E_grid, dat_name2.c_str());

    for(unsigned int r=0;r<layer.r_dim;++r){
      for(unsigned int z=0;z<layer.z_dim;++z){
        diff_fields[r][z]=fabs(efield[r][z]-std_efield[r][z]);}}

    /******** DIFFERENCE IN EFIELD OUTPUT ************/
    FILE  *file3 = fopen(dat_name3.c_str(), "w");
    if ( !file3 ) 
    { printf("could not open %s\n",dat_name3.c_str()); 
      exit(EXIT_FAILURE); }
    for ( unsigned int r = 0; r < layer.r_dim; ++r){ 
      for (unsigned int z = 0; z < layer.z_dim; ++z){
        fprintf(file3, " % .5e", diff_fields[r][z]);
      }
      fprintf(file3, "\n");
    }
    fclose(file3);
  }
#endif
}

std::vector<std::vector<double>> out_std_efield ( Field E_grid[], std::string dat_name2 ) {
  FILE  *file2 = fopen(dat_name2.c_str(), "w");
  if ( !file2 )
  { printf("could not open %s\n",dat_name2.c_str()); 
    exit(EXIT_FAILURE); }

  GridLayer& layer=global_grid.layers[ELECTRONS];
  unsigned int rmax, zmax;
  rmax = layer.r_dim-1;
  zmax = layer.z_dim-1;
  int NZ = global_grid.mesh_z_dim;
  double f = 1./1.;

  /*************** normal ****************************/
  std::vector<std::vector<double>> stdfield        = {{0.0}};
  /************ cell face field arrays ***************/
  std::vector<std::vector<double>> left_stdfield   = {{0.0}};
  std::vector<std::vector<double>> right_stdfield  = {{0.0}};
  std::vector<std::vector<double>> top_stdfield    = {{0.0}};
  std::vector<std::vector<double>> bottom_stdfield = {{0.0}};

  /*************** normal ****************************/
  stdfield.resize(layer.r_dim);
  /************ cell face field arrays ***************/
  left_stdfield.resize(layer.r_dim);
  top_stdfield.resize(layer.r_dim);
  right_stdfield.resize(layer.r_dim);
  bottom_stdfield.resize(layer.r_dim);

  for (unsigned int r = 0; r <= rmax; ++r) {

    /*************** normal **************************/
    stdfield[r].resize(layer.z_dim);
    /************ cell face field arrays *************/
    top_stdfield[r].resize(layer.z_dim);
    bottom_stdfield[r].resize(layer.z_dim);
    left_stdfield[r].resize(layer.z_dim);
    right_stdfield[r].resize(layer.z_dim);

    for (unsigned int z = 0; z <= zmax; ++z) {

      /************ cell fac adaptation **************/
      // bottom left grid point (main)
      bottom_stdfield[r][z]  =  (1./2.)*E_grid[ r   *NZ+ z   ].r*f;
      left_stdfield[r][z]    =  (1./2.)*E_grid[ r   *NZ+ z   ].z*f;
      // top right grid point - corner
      top_stdfield[r][z]     =  (1./2.)*E_grid[(r+1)*NZ+(z+1)].r*f;
      right_stdfield[r][z]   =  (1./2.)*E_grid[(r+1)*NZ+(z+1)].z*f;
      // top left grid point - corner
      left_stdfield[r][z]   +=  (1./2.)*E_grid[(r+1)*NZ+ z   ].z*f;
      top_stdfield[r][z]    +=  (1./2.)*E_grid[(r+1)*NZ+ z   ].r*f;
      // bottom right grid point - corner
      right_stdfield[r][z]  +=  (1./2.)*E_grid[ r   *NZ+(z+1)].z*f;
      bottom_stdfield[r][z] +=  (1./2.)*E_grid[ r   *NZ+(z+1)].r*f;

      /************ summarization ********************/
      stdfield[r][z] = ( right_stdfield[r][z] + left_stdfield[r][z] 
                        + top_stdfield[r][z] + bottom_stdfield[r][z] );

      iiiprintf( ">> matrix2fem efield: r=%i/%i, z=%i/%i; "
                 "stdfield[r][z]=%g, top_stdfield=%g,"
                 " bottom_stdfield=%g, left_stdfield=%g, right_stdfield=%g; "
                 "topleftgrid.r=%g, toprightgrid.r=%g, bottomleftgrid.r=%g, "
                 "bottomrightgrid.r=%g, "
                 "topleftgrid.z=%g, bottomleftgrid.z=%g, bottomrightgrid.z=%g, toprightgrid.z=%g\n",
                 r, rmax, z, zmax, stdfield[r][z], top_stdfield[r][z], bottom_stdfield[r][z],
                 left_stdfield[r][z], right_stdfield[r][z],
                 E_grid[(r+1)*NZ+z].r, E_grid[(r+1)*NZ+(z+1)].r,
                 E_grid[r*NZ+z].r, E_grid[r*NZ+(z+1)].r,
                 E_grid[(r+1)*NZ+z].z, E_grid[r*NZ+z].z,
                 E_grid[r*NZ+(z+1)].z, E_grid[(r+1)*NZ+(z+1)].z );
      
      fprintf ( file2, " % .5e", stdfield[r][z] );
    }
    fprintf ( file2, "\n" );
  }
  fclose ( file2 );
  return stdfield;
}

void out_vector_femfield ( std::string dat_name ) {

}

void out_vector_stdfield ( Field E_grid[], std::string dat_name ) {
#if USE_FEM_SOLVER
  int mpi_rank=get_rank();
  if (mpi_rank==0) {
    FILE  *file = fopen(dat_name.c_str(), "w");
    if ( !file )
    { printf("could not open %s\n",dat_name.c_str()); i
      exit(EXIT_FAILURE); }
    
    GridLayer& layer=global_grid.layers[ELECTRONS];
    unsigned int rmax, zmax;
    rmax = layer.r_dim-1;
    zmax = layer.z_dim-1;
    int NZ = global_grid.mesh_z_dim;
    double f = 1./1.;
    double max_norm = 0.0;
    
    /*************** normal ****************************/
    std::vector<std::vector<double>> f_stdfield = {{0.0}};
    /************ cell face field arrays ***************/
    std::vector<std::vector<double>> z_stdfield = {{0.0}};
    std::vector<std::vector<double>> r_stdfield = {{0.0}};
    
    /*************** normal ****************************/
    f_stdfield.resize(layer.r_dim);
    /************ cell face field arrays ***************/
    z_stdfield.resize(layer.r_dim);
    r_stdfield.resize(layer.r_dim);
    
    for (unsigned int r = 0; r <= rmax; ++r) {
    
      /*************** normal **************************/
      f_stdfield[r].resize(layer.z_dim);
      /************ cell face field arrays *************/
      z_stdfield[r].resize(layer.z_dim);
      r_stdfield[r].resize(layer.z_dim);
    
      for (unsigned int z = 0; z <= zmax; ++z) {
        
        /************ cell fac adaptation **************/
        // bottom left grid point (main)
        r_stdfield[r][z] =  (1./2.)*E_grid[ r   *NZ+ z   ].r*f;
        z_stdfield[r][z] =  (1./2.)*E_grid[ r   *NZ+ z   ].z*f;
        // top right grid point - corner
        r_stdfield[r][z] += (1./2.)*E_grid[(r+1)*NZ+(z+1)].r*f;
        z_stdfield[r][z] += (1./2.)*E_grid[(r+1)*NZ+(z+1)].z*f;
        // top left grid point - corner
        z_stdfield[r][z] += (1./2.)*E_grid[(r+1)*NZ+ z   ].z*f;
        r_stdfield[r][z] += (1./2.)*E_grid[(r+1)*NZ+ z   ].r*f;
        // bottom right grid point - corner
        z_stdfield[r][z] += (1./2.)*E_grid[ r   *NZ+(z+1)].z*f;
        r_stdfield[r][z] += (1./2.)*E_grid[ r   *NZ+(z+1)].r*f;
    
        /************ summarization ********************/
        f_stdfield[r][z] = sqrt(SQU(r_stdfield[r][z]) +
                                SQU(z_stdfield[r][z]));
      
        /*********** need maximum to calculate norm ****/
        /*********** for arrows length/color ***********/
        max_norm = std::max( max_norm, f_stdfield[r][z] );
        
      }
    }
    
    for (unsigned int r = 0; r <= rmax; ++r) {
      for (unsigned int z = 0; z <= zmax; ++z) {
        fprintf ( file, " % .5e % .5e % .5e % .5e % .5e ",
                  z, r,                                                 // x&y
                  z_stdfield[r][z]/max_norm, r_stdfield[r][z]/max_norm, // dx&dy
                  f_stdfield[r][z] );                                   // color, e.g. strength
      } fprintf ( file, "\n" );
    }   fclose ( file );

  }
#endif
}
