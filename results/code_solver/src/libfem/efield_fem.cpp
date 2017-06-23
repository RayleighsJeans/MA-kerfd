#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string>
#include <iostream>
#include <fstream>

#include "grid.h"
#include "pic.h"
#include "fem_solver.h"
#include "debug_printing.h"
#define DEBUG_LEVEL DEBUG_INFO_3

void efield_matrix2fem ( Field E_grid[] ) {
#if USE_FEM_SOLVER

  // deposit electric field of matrix solver method 
  // to FEM cell properties according to prerequisities 
  // for solvers time evolution

  // necessary cell grab
  GridLayer& layer=global_grid.layers[ELECTRONS];
  unsigned int NZ, rmax, zmax;
  NZ = global_grid.mesh_z_dim;
  rmax = layer.r_dim-1;
  zmax = layer.z_dim-1;
  double f = 1./1.;

  for ( unsigned int r = 0; r <= rmax; ++r){
    for (unsigned int z = 0; z <= zmax; ++z){
      Cell& cell = layer.get_cell(r,z);
     
      // for each cell face to grid points have to be accounted for
      // choose for coefficients of correct weighting for electrid field of grid to cell
      // each cell neighbouring a grid point has to recieve equal portions of electric field
      // so no ghost field is calculated and neither is lost

      // bottom left grid point (main)
      cell.fem.efield_bottom = (1./2.)*E_grid[ r   *NZ+ z   ].r*f;
      cell.fem.efield_left   = (1./2.)*E_grid[ r   *NZ+ z   ].z*f;
      // top right grid point - corner
      cell.fem.efield_top    = (1./2.)*E_grid[(r+1)*NZ+(z+1)].r*f;
      cell.fem.efield_right  = (1./2.)*E_grid[(r+1)*NZ+(z+1)].z*f;
      // top left grid point - corner
      cell.fem.efield_left   += (1./2.)*E_grid[(r+1)*NZ+ z   ].z*f;
      cell.fem.efield_top    += (1./2.)*E_grid[(r+1)*NZ+ z   ].r*f;
      // bottom right grid point - corner
      cell.fem.efield_right  += (1./2.)*E_grid[ r   *NZ+(z+1)].z*f;
      cell.fem.efield_bottom += (1./2.)*E_grid[ r   *NZ+(z+1)].r*f;

      iiiprintf( ">> matrix2fem efield: r=%i/%i, z=%i/%i; "
                 "topfield=%g, bottomfield=%g, leftfield=%g, rightrield=%g; "
                 "topleftgrid.r=%g, toprightgrid.r=%g, bottomleftgrid.r=%g, bottomrightgrid.r=%g, "
                 "topleftgrid.z=%g, bottomleftgrid.z=%g, bottomrightgrid.z=%g, toprightgrid.z=%g\n",
                 r, rmax, z, zmax,
                 cell.fem.efield_top, cell.fem.efield_bottom, cell.fem.efield_left, cell.fem.efield_right,
                 E_grid[(r+1)*NZ+z].r, E_grid[(r+1)*NZ+(z+1)].r, E_grid[r*NZ+z].r, E_grid[r*NZ+(z+1)].r,
                 E_grid[(r+1)*NZ+z].z, E_grid[r*NZ+z].z, E_grid[r*NZ+(z+1)].z, E_grid[(r+1)*NZ+(z+1)].z );
    }
  }

  out_efield( "out/fem_efield_relax.dat",
              "out/stand_efield_relax.dat",
              "out/diff_fields_relax.dat",
              E_grid,
              "out/old_field_relax.dat" );

#endif
}

void efield_fem ( ) {
#if USE_FEM_SOLVER
  /********************************************************************************
  // !!! ATTENTION ATTENTION ATTENTION ATTENTION ATTENTION ATTENTION !!!
  // EFIELD WILL BE SUMMED UP FOR ALL CELL FACES, SO EACH CELL FACE
  // ACCUMULATES THE SAME FIELD, WHICH TRANSLATES TO THE NEEDLESS SUMMERIZATION
  // INSIDE THE PUSH ODER THE OLD FIELD DUMP PROCESS
  // EASES THE TRANSFORMATION OF THE MATRIX/GRID FIELD OF SLU/SOR TO CELL PROPERTIES
  **********************************************************************************/  

  double boundr, boundz;
         boundr = boundz = 0.0;
  double Ua_ext, Ua_SB, qe, fr, fz;
         Ua_ext = Ua_SB = qe = fr = fz = 0.0;
         qe = SQU(dt)/(2.0*Ncell1);
  double scale_pot = 0.5*SQU(dt/dr)/(T_e0);
  double scale_efield = 0.5*dr_0*scale_pot;
         Ua_SB = -200*0.5*SQU(dt/dr)/T_e0;
  iiiprintf(">> efield calculation: self bias"
            " voltage used is %g, scale pot = %g,"
            " scale_efield =%g and qe=%g\n",
            Ua_SB, scale_pot, scale_efield, qe);

  // necessary cell grab
  GridLayer& layer=global_grid.layers[ELECTRONS];
  for ( unsigned int r = 0; r < layer.r_dim; ++r){
    for (unsigned int z = 0; z < layer.z_dim; ++z){
      Cell& cell = layer.get_cell(r,z);

      // scaling factors for ok measures
      if (r==0) { fz = -6.0*qe; fr = -qe; } else if ( r == layer.r_dim-1 ) {
      fr = -6.0*qe/(3.0*layer.r_dim-1.0); fz = -qe/(layer.r_dim-1); } else {
      fz = -qe/r; fr = -qe/(r+1); }
      iiiprintf(">> efield calculation: r/z=%i/%i, "
                "fr = %g, fz = %g scaling\n",
                r, z, fr, fz);
                
      if ( Eps[r][z] == -100 ) {

#if GEOM_TYPE_RF
        iiiprintf(">> efield calculation: "
                  "external/const pot at r %i & z %i\n", r, z);
        if ( phi_bound[r*(global_grid.z_dim+1)+z] != -1 ) // accounting for cathode/anode
        { Ua_ext = phi_bound[r*(global_grid.z_dim+1)+z] * 
                   sin(2*PI*(double)nstep/(f_RF))+Ua_SB;
          iiiprintf(">> efield calculation: external pot"
                    " at r=%i,z=%i; Ua_ext=%g, phi_bound=%f\n",
                    r, z, Ua_ext, phi_bound[r*(global_grid.z_dim+1)+z] ); }
        
        if ( r < layer.r_dim-1 ) {                        // differential quotient for cathode/anode

          iiiprintf(">> efield calculation: differrential quotient"
                    " for cath/anode at r=%i/%i, z=%i/%i\n",
                    r, layer.r_dim-1, z, layer.z_dim-1);
        } else {                                          // differential quotient for top boundary
          iiiprintf(">> efield calculation: differrential quotient"
                    " for top at r=%i/%i, z=%i/%i\n",
          r, layer.r_dim-1, z, layer.z_dim-1);
        }
#endif
        iiiprintf(">> efield calculation: break from loop at r=%i/%i, z=%i/%i, Eps[r][z]=%i\n",
        r, layer.r_dim-1, z, layer.z_dim-1, Eps[r][z]);
     
       } else if ( (Eps[r][z] == -300) || (Eps[r][z] == -200) ) { 
      
        boundr = fabs(((double)layer.r_dim-1)-(double)r);
        boundz = fabs(((double)layer.z_dim-1)-(double)z);
        iiiprintf(">> efield calculation: const "
                  "efield 0 at r=%i, z=%i, boundr=%g, boundz=%g \n",
                  r, z, boundr, boundz );
        // accouting for const efield at walls/top
        // left or right of domain; z component is zero
        if ( boundr >= boundz ) {
          if ( z == 0 ) { 
            cell.fem.efield_left   = 0.0;
          } else if ( z == layer.z_dim-1 )  {
            cell.fem.efield_right  = 0.0;
          } else {
            cell.fem.efield_right =
            cell.fem.efield_left = 0.0;
          }
        // top or bottom of domain; r component is zero
        } else {
          if ( r == 0 ) { 
            cell.fem.efield_bottom   = 0.0;
          } else if ( r == layer.r_dim-1 )  {
            cell.fem.efield_top  = 0.0;
          } else {
            cell.fem.efield_top =
            cell.fem.efield_bottom = 0.0;
          }
        }
        iiiprintf(">> efield calculation: break from loop at r=%i/%i, z=%i/%i, Eps[r][z]=%i\n",
        r, layer.r_dim-1, z, layer.z_dim-1, Eps[r][z]);
      
      } else {

        cell.fem.efield_top    =  (- ( 2.*TWOPI*dt*cell.fem.current_top ) +
                                  cell.fem.oldfield_top );
        cell.fem.efield_bottom =  (- ( 2.*TWOPI*dt*cell.fem.current_bottom ) +
                                  cell.fem.oldfield_bottom );
        cell.fem.efield_left   =  (- ( 2.*TWOPI*dt*cell.fem.current_left ) +
                                  cell.fem.oldfield_left );
        cell.fem.efield_right  =  (- ( 2.*TWOPI*dt*cell.fem.current_right ) +
                                  cell.fem.oldfield_right );

        if ( z > 0 )
        { Cell& leftcell = layer.get_cell(r,z-1);
          cell.fem.efield_left += - ( 2.*TWOPI*dt*cell.fem.current_right );
          iiiprintf(">> efield calculation : leftcell found, efield_left=%g, r=%d/%d, z=%d/%d\n",
                     cell.fem.efield_left, r, layer.r_dim-1, z-1, layer.z_dim-1); } 
        if ( r > 0 )
        { Cell& bottomcell = layer.get_cell(r-1,z);
          cell.fem.efield_bottom += - ( 2.*TWOPI*dt*bottomcell.fem.current_top );
          iiiprintf(">> efield calculation : bottomcell found, efield_bottom=%g, r=%d/%d, z=%d/%d\n",
                    cell.fem.efield_bottom, r-1, layer.r_dim-1, z, layer.z_dim-1); }
        if ( r < layer.r_dim-1 )
        { Cell& topcell = layer.get_cell(r+1,z);
          cell.fem.efield_top += - ( 2.*TWOPI*dt*topcell.fem.current_bottom );
          iiiprintf(">> efield calculation : topcell found, efield_top=%g, r=%d/%d, z=%d/%d\n",
                    cell.fem.efield_top, r+1, layer.r_dim-1, z, layer.z_dim-1); }
        if ( z < layer.z_dim-1 )
        { Cell& rightcell = layer.get_cell(r,z+1);
          cell.fem.efield_right += - ( 2.*TWOPI*dt*rightcell.fem.current_left );
          iiiprintf(">> efield calculation : rightcell found, efield_right=%g, r=%d/%d, z=%d/%d\n",
                    cell.fem.efield_right, r, layer.r_dim-1, z+1, layer.z_dim-1); }
        iiiprintf(">> efield calculation : r=%i/%i, z=%i/%i; fields top=%g, bottom=%g, left=%g, right=%g; "
                  " currents top=%g, bottom=%g, left=%g, right=%g\n",
                  r, layer.r_dim-1, z, layer.z_dim-1, cell.fem.efield_top, cell.fem.efield_bottom,
                  cell.fem.efield_left, cell.fem.efield_right, cell.fem.current_top,
                  cell.fem.current_bottom, cell.fem.current_left, cell.fem.current_right ); 
      }
    }
  }
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
  iiprintf(">> reset fem cell diagnostic "
           "arrays ... successfull!\n");
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
  iiiprintf(">> reset fem old field and deposition"
            " of new field as old ... successfull!\n");
#endif
}
