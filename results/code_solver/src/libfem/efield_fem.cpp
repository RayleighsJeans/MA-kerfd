#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <string>
#include <iostream>
#include <fstream>

#include "pic.h"
#include "grid.h"

#include "fem_solver.h"
#include "debug_printing.h"
#include "fem_debug.h"

void efield_fem ( ) {
#if USE_FEM_SOLVER

  // necessary cell grab
  GridLayer& layer=global_grid.layers[ELECTRONS];
  for ( unsigned int r = 0; r < layer.r_dim; ++r){
    for (unsigned int z = 0; z < layer.z_dim; ++z){
      Cell& cell = layer.get_cell(r,z);

      cell.fem.efield_top    =  - ( dt*cell.fem.current_top ) +
                                                     cell.fem.oldfield_top;
      cell.fem.efield_bottom =  - ( dt*cell.fem.current_bottom ) +
                                                     cell.fem.oldfield_bottom;
      cell.fem.efield_left   =  - ( dt*cell.fem.current_left ) +
                                                     cell.fem.oldfield_left;
      cell.fem.efield_right  =  - ( dt*cell.fem.current_right ) +
                                                     cell.fem.oldfield_right;


      if ( z > 0 ) {
        Cell& leftcell = layer.get_cell(r,z-1);
        cell.fem.efield_left  +=  - ( dt*leftcell.fem.current_right ) +
                                                       leftcell.fem.oldfield_right;
      } else if ( r > 0 ) {
        Cell& bottomcell = layer.get_cell(r-1,z);
        cell.fem.efield_bottom+=  - ( dt*bottomcell.fem.current_top ) +
                                                       bottomcell.fem.oldfield_top;
      } else if ( r < layer.r_dim-1 ) {
        Cell& topcell = layer.get_cell(r+1,z);
        cell.fem.efield_top   +=  - ( dt*topcell.fem.current_bottom ) +
                                                       topcell.fem.oldfield_bottom;
      } else if ( z < layer.z_dim-1 ) {
        Cell& rightcell = layer.get_cell(r,z+1);
        cell.fem.efield_right +=  - ( dt*rightcell.fem.current_left ) +
                                                       rightcell.fem.oldfield_left;
      }

    }
  }
  store_old_field ( );

#endif
}

void store_old_field ( ) {
#if USE_FEM_SOLVER

  // necessary cell grab
  GridLayer& layer=global_grid.layers[ELECTRONS];
  for ( unsigned int r = 0; r < layer.r_dim; ++r){
    for (unsigned int z = 0; z < layer.z_dim; ++z){
      Cell& cell = layer.get_cell(r,z);

      // reset oldfield
       cell.clear_oldfield();

      // store current field as new field
      cell.fem.oldfield_top    = cell.fem.efield_top;
      cell.fem.oldfield_bottom = cell.fem.efield_bottom;
      cell.fem.oldfield_left   = cell.fem.efield_left;
      cell.fem.oldfield_right  = cell.fem.efield_right;
    }
  }
  
  iprintf(">> reset fem old field and deposition of new field as old ... successfull!");

#endif
}

void reset_fem ( ) {
#if USE_FEM_SOLVER

  // necessary cell grab
  GridLayer& layer=global_grid.layers[ELECTRONS];
  for ( unsigned int r = 0; r < layer.r_dim; ++r){
    for (unsigned int z = 0; z < layer.z_dim; ++z){
      Cell& cell = layer.get_cell(r,z);
  
      // reset fem diagnostic arrays
      cell.clear_fem();
        
    }
  }

  iprintf(">> reset fem cell diagnostic arrays ... successfull!");

#endif
}
