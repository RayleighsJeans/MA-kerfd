#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string>
#include <iostream>
#include <fstream>

#include "pic.h"
#include "fem_solver.h"
#include "debug_printing.h"
#include "fem_debug.h"

void efield_fem ( ) {
#if USE_FEM_SOLVER
  
  double Ua_ext, Ua_SB;
  double scale_pot = 0.5*SQU(dt/dr)/(T_e0);
         Ua_SB = -200*0.5*SQU(dt/dr)/T_e0;
  iiiprintf(">> efield calculation: self bias"
            " voltage used is %g, scale pot = %g\n",
            Ua_SB, scale_pot);
  
  // necessary cell grab
  GridLayer& layer=global_grid.layers[ELECTRONS];
  for ( unsigned int r = 0; r < layer.r_dim; ++r){
    for (unsigned int z = 0; z < layer.z_dim; ++z){
      Cell& cell = layer.get_cell(r,z);
      cell.fem.efield_top    =  (- ( dt*cell.fem.current_top ) +
                                cell.fem.oldfield_top);
      cell.fem.efield_bottom =  (- ( dt*cell.fem.current_bottom ) +
                                cell.fem.oldfield_bottom);
      cell.fem.efield_left   =  (- ( dt*cell.fem.current_left ) +
                                cell.fem.oldfield_left);
      cell.fem.efield_right  =  (- ( dt*cell.fem.current_right ) +
                                cell.fem.oldfield_right);

      if ( Eps[r][z] == -100 ) {
        // ACCOUNTING FOR APPLIED VOLTAGE
        if( phi_bound[r*(global_grid.z_dim+1)+z] != -1 ){
          Ua_ext = phi_bound[r*(global_grid.z_dim+1)+z] * 
                   sin(2*PI*(double)nstep/(f_RF))+Ua_SB;
          iiiprintf(">> efield calculation: external pot at r=%i,z=%i;i Ua_ext=%g, phi_bound=%f\n",
                    r, z, Ua_ext, phi_bound[r*(global_grid.z_dim+1)+z] );
        }
        cell.fem.efield_left = 0.5*( (cell.fem.area_weighted_charge*scale_pot) -
                               Ua_ext );
        iiiprintf(">> efield calculation: additional field at r=%i,z=%i; efield_left=%g, Eps[r][z]=%g\n",
                  r, z, cell.fem.efield_left, Eps[r][z] );
      } else if ( (Eps[r][z] == -300) || (Eps[r][z] == -200) ) {
        // ACCOUNTING FOR CONST EFIELD AT WALLS
        cell.fem.efield_right = 0.0;
        iiiprintf(">> efield calculation: const efield 0 at r=%i, z=%i\n",
                  r, z );
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
  iiprintf(">> reset fem old field and deposition of new field as old ... successfull!\n");

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
  iiprintf(">> reset fem cell diagnostic arrays ... successfull!\n");

#endif
}
