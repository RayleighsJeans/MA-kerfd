#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <string.h>
#include <math.h>
#include <cmath>
#include <vector>
#include <typeinfo>

#include <string>
#include <iostream>
#include <fstream>

#include "grid.h"
#include "pic.h"
#include "boundary_check.h"
#include "fem_solver.h"
#include "debug_printing.h"
#define DEBUG_LEVEL DEBUG_INFO_3



void efield_matrix2fem ( Field E_grid[] ) {
#if USE_FEM_SOLVER

  /*************************************************************/
  // deposit electric field of matrix solver method 
  // to FEM cell properties according to prerequisities 
  // for solvers time evolution
  // for each cell face to grid points
  // have to be accounted for
  // choose for coefficients of correct weighting for
  // electrid field of grid to cell. each cell neighbouring
  //a grid point has to recieve equal portions of electric field
  // so no ghost field is calculated and neither is lost
  /*************************************************************/

  GridLayer& layer=global_grid.layers[ELECTRONS];
  unsigned int NZ, rmax, zmax;
  NZ = global_grid.mesh_z_dim;
  rmax = layer.r_dim-1;
  zmax = layer.z_dim-1;
  double f = 1./1.;

  for ( unsigned int r = 0; r <= rmax; ++r){
    for (unsigned int z = 0; z <= zmax; ++z){
      Cell& cell = layer.get_cell(r,z);
     
      // bottom left grid point (main)
      cell.fem.efield_bottom  = (1./2.)*E_grid[ r   *NZ+ z   ].r*f;
      cell.fem.efield_left    = (1./2.)*E_grid[ r   *NZ+ z   ].z*f;
      // top right grid point - corner
      cell.fem.efield_top     = (1./2.)*E_grid[(r+1)*NZ+(z+1)].r*f;
      cell.fem.efield_right   = (1./2.)*E_grid[(r+1)*NZ+(z+1)].z*f;
      // top left grid point - corner
      cell.fem.efield_left   += (1./2.)*E_grid[(r+1)*NZ+ z   ].z*f;
      cell.fem.efield_top    += (1./2.)*E_grid[(r+1)*NZ+ z   ].r*f;
      // bottom right grid point - corner
      cell.fem.efield_right  += (1./2.)*E_grid[ r   *NZ+(z+1)].z*f;
      cell.fem.efield_bottom += (1./2.)*E_grid[ r   *NZ+(z+1)].r*f;

      iiiprintf( ">> matrix2fem efield: r=%i/%i, z=%i/%i; "
                 "topfield=%g, bottomfield=%g, leftfield=%g, "
                 "rightrield=%g; topleftgrid.r=%g, toprightgrid.r=%g,"
                 " bottomleftgrid.r=%g, bottomrightgrid.r=%g, "
                 "topleftgrid.z=%g, bottomleftgrid.z=%g,"
                 " bottomrightgrid.z=%g, toprightgrid.z=%g\n",
                 r, rmax, z, zmax,
                 cell.fem.efield_top, cell.fem.efield_bottom,
                 cell.fem.efield_left, cell.fem.efield_right,
                 E_grid[(r+1)*NZ+z].r, E_grid[(r+1)*NZ+(z+1)].r,
                 E_grid[r*NZ+z].r, E_grid[r*NZ+(z+1)].r,
                 E_grid[(r+1)*NZ+z].z, E_grid[r*NZ+z].z,
                 E_grid[r*NZ+(z+1)].z, E_grid[(r+1)*NZ+(z+1)].z );
    } // z -> layer.z_dim-1
  } // r -> layer.r_dim-1

  out_efield ( "out/fem_efield_relax.dat",
               "out/stand_efield_relax.dat",
               "out/diff_fields_relax.dat",
               E_grid, "out/old_field_relax.dat" );
  out_vector_femfield ( "out/vector_femfield_relax.dat" );
  out_vector_stdfield ( E_grid, "out/vector_stdfield_relax.dat" );

#endif
}

void fem_eps2cell ( ) {
#if USE_FEM_SOLVER

  GridLayer& layer = global_grid.layers[ELECTRONS];
  unsigned int rmax, zmax;
             rmax = layer.r_dim-1;
             zmax = layer.z_dim-1;
  
  for ( unsigned int r = 0; r < layer.r_dim; ++r){
    for (unsigned int z = 0; z < layer.z_dim; ++z){
      Cell& cell = layer.get_cell(r,z);
      cell.fem.epsilon = eps_grid2face ( r, z, rmax, zmax );
      
      iiiprintf ( ">> efield calculation: r/z=%i/%i, "
                  "epsvec(0-3)=%i,%i,%i,%i while "
                  " Eps(gridmatrix)=%g,%g,%g,%g \n",
                  r, z, cell.fem.epsilon[0], cell.fem.epsilon[1],
                  cell.fem.epsilon[2],cell.fem.epsilon[3],
                  Eps[r][z], Eps[r+1][z],
                  Eps[r+1][z+1], Eps[r][z+1] );

    }
  }
  
  iiiprintf ( ">> eps2cell: moved Epsilon"
              " from gridpoints to cell properties! \n" );

#endif
}

// ernumeration of cell faces for eps grid 
// 2 cell wall algorithm
//
//(r+1,z)                (r+1,z+1)
//topleft                 topright  
//  **************************
//  *           1            *   
//  *                        * 
//  *                        * 
//  *                        * 
//  * 0                    2 *
//  *                        *
//  *                        *
//  *                        *
//  *           3            *
//  **************************
//bottomleft             bottomright
// (r,z)                    (r,z+1)
//

std::vector<int> eps_grid2face ( unsigned int r, unsigned int z,
                                 unsigned int rmax, unsigned int zmax ) {
#if USE_FEM_SOLVER

  std::vector<int> face_eps;
    face_eps.resize(4); face_eps = { 1 };

  face_eps[0] = std::round(Eps[r+1][z]);
  face_eps[1] = std::round(Eps[r+1][z+1]);
  face_eps[2] = std::round(Eps[r+1][z+1]);
  face_eps[3] = std::round(Eps[r][z+1]);

  if ( r == rmax ) {
    face_eps[2] = std::round(Eps[r][z+1]);
    face_eps[0] = std::round(Eps[r][z]);
  } else if ( r == 0 ) {
    face_eps[2] = std::round(Eps[r+1][z+1]);
    face_eps[0] = std::round(Eps[r+1][z]);
  }

  if ( z == zmax ) {
    face_eps[1] = std::round(Eps[r+1][z]);
    face_eps[3] = std::round(Eps[r][z]);
  } else if ( z == 0 ) {
    face_eps[1] = std::round(Eps[r+1][z+1]);
    face_eps[3] = std::round(Eps[r][z+1]);
  }

  if (((face_eps[1] == -100) && (r != rmax)) || 
       (face_eps[1] == -300)                 || 
       (face_eps[1] == -200)                 ) {
     face_eps[1] = 1; }
  
  if (((face_eps[3] == -100) && (r != 0   )) || 
       (face_eps[3] == -200)                 || 
       (face_eps[3] == -300)                 || 
       (face_eps[3] == -400)                 ) {
    face_eps[3] = 1; }

  if (face_eps[0] == -400) {
    face_eps[0] = 1; }

  if (face_eps[2] == -400) {
    face_eps[2] = 1; }

  return face_eps;
  
#endif
}

void efield_constpot ( unsigned int r, unsigned int z,
                       std::vector<int> &epsvec, double Ua_SB,
                       Cell& cell, BoundarySegment& polygon,
                       std::vector<unsigned int> &queue ) {
#if USE_FEM_SOLVER
  double Ua_ext = 0.0;

  // searching for boundary condition to const pot
  if ( epsvec[0] == -100 ) {
    queue[0] = 0;
#if GEOM_TYPE_RF
    if ( phi_bound[r*(global_grid.z_dim+1)+z] != -1 )
    { Ua_ext = phi_bound[r*(global_grid.z_dim+1)+z] * 
               sin(2*PI*(double)nstep/(f_RF))+Ua_SB; }
    else if ( phi_bound[(r+1)*(global_grid.z_dim+1)+z] != -1 )
    { Ua_ext = phi_bound[(r+1)*(global_grid.z_dim+1)+z] * 
               sin(2*PI*(double)nstep/(f_RF))+Ua_SB; }
#endif
    iiiprintf ( ">> efield calculation: external pot left"
                " at r=%i,z=%i; Ua_ext=%g, phi_bound_down=%g,"
                " phi_bound_top=%g,"
                " differential quotient for electrodes/top\n",
                r, z, Ua_ext, phi_bound[r*(global_grid.z_dim+1)+z],
                phi_bound[(r+1)*(global_grid.z_dim+1)+z] );
  }  // in epsvec[0] -100
 
  if ( epsvec[1] == -100 ) {
    queue[1] = 0;
#if GEOM_TYPE_RF
    if ( phi_bound[(r+1)*(global_grid.z_dim+1)+z] != -1 )
    { Ua_ext = phi_bound[(r+1)*(global_grid.z_dim+1)+z] * 
              sin(2*PI*(double)nstep/(f_RF))+Ua_SB; }
    else if ( phi_bound[(r+1)*(global_grid.z_dim+1)+z+1] != -1 )
    { Ua_ext = phi_bound[(r+1)*(global_grid.z_dim+1)+z+1] * 
              sin(2*PI*(double)nstep/(f_RF))+Ua_SB; }
#endif
    iiiprintf ( ">> efield calculation: external pot top"
                " at r=%i,z=%i; Ua_ext=%g, phi_bound_left=%g,"
                " phi_bound_right=%g,"
                " differential quotient for electrodes/top\n",
                r, z, Ua_ext, phi_bound[(r+1)*(global_grid.z_dim+1)+z],
                phi_bound[(r+1)*(global_grid.z_dim+1)+z+1] );
  }  // in epsvec[1] -100

  if ( epsvec[2] == -100 ) {
    queue[2] = 0;
#if GEOM_TYPE_RF
    if ( phi_bound[(r+1)*(global_grid.z_dim+1)+z+1] != -1 )
    { Ua_ext = phi_bound[(r+1)*(global_grid.z_dim+1)+(z+1)] * 
              sin(2*PI*(double)nstep/(f_RF))+Ua_SB; }
    else if ( phi_bound[r*(global_grid.z_dim+1)+z+1] != -1 )
    { Ua_ext = phi_bound[r*(global_grid.z_dim+1)+z+1] * 
              sin(2*PI*(double)nstep/(f_RF))+Ua_SB; }
#endif
    iiiprintf ( ">> efield calculation: external pot right"
                " at r=%i,z=%i; Ua_ext=%g, phi_bound_top=%g,"
                " phi_bound_bottom=%g,"
                " differential quotient for electrodes/top\n",
                r, z, Ua_ext, phi_bound[(r+1)*(global_grid.z_dim+1)+z+1],
                phi_bound[r*(global_grid.z_dim+1)+z+1] );
  }  // in epsvec[2] -100

  if ( epsvec[3] == -100 ) {
    queue[3] = 0;
#if GEOM_TYPE_RF
    if ( phi_bound[r*(global_grid.z_dim+1)+z] != -1 )
    { Ua_ext = phi_bound[r*(global_grid.z_dim+1)+z] * 
              sin(2*PI*(double)nstep/(f_RF))+Ua_SB; }
    else if ( phi_bound[r*(global_grid.z_dim+1)+z+1] != -1 )
    { Ua_ext = phi_bound[r*(global_grid.z_dim+1)+z+1] * 
              sin(2*PI*(double)nstep/(f_RF))+Ua_SB; }
#endif
    iiiprintf ( ">> efield calculation: external pot bottom"
                " at r=%i,z=%i; Ua_ext=%g, phi_bound_left=%g,"
                " phi_bound_right=%g,"
                " differential quotient for electrodes/top\n",
                r, z, Ua_ext, phi_bound[r*(global_grid.z_dim+1)+z],
                phi_bound[r*(global_grid.z_dim+1)+z+1] );
  }  // in epsvec[3] -100

#endif
}

void efield_constfield ( unsigned int r, unsigned int z,
                         std::vector<int> &epsvec,
                         Cell& cell, BoundarySegment& polygon,
                         std::vector<unsigned int> &queue ) {
#if USE_FEM_SOLVER

  if ( epsvec[0] == -300 ) {
    iiiprintf ( ">> efield calculation: const "
                "efield left %g at r=%i, z=%i, epsvec[1]=%i\n",
                polygon.efield[0], r, z, epsvec[0] );
    cell.fem.efield_left = polygon.efield[0]; queue[0] = 0; }

  if ( epsvec[1] == -400 ) {
    iiiprintf ( ">> efield calculation: const "
                "efield top %g at r=%i, z=%i, epsvec[1]=%i\n",
                polygon.efield[1], r, z, epsvec[1] );
    cell.fem.efield_top = polygon.efield[1]; queue[1] = 0; }

  if ( epsvec[2] == -200 ) {
    iiiprintf ( ">> efield calculation: const "
                "efield right %g at r=%i, z=%i, epsvec[2]=%i\n",
                polygon.efield[2], r, z, epsvec[2] );
    cell.fem.efield_right = polygon.efield[2]; queue[2] = 0; }
  
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

        GridLayer& layer = global_grid.layers[ELECTRONS];
      unsigned int rmax, zmax;
                   rmax = layer.r_dim-1;
                   zmax = layer.z_dim-1;
std::vector<unsigned int>
                   queue;
                   queue.resize(4); queue = { 1 };
            double Ua_SB, qe, fr, fz;
                   Ua_SB = qe = fr = fz = 0.0;
                   qe = SQU(dt)/(2.0*Ncell1);
            double scale_pot = 0.5*SQU(dt/dr)/(T_e0);
            double scale_efield = 0.5*dr_0*scale_pot;
                   Ua_SB = -200*0.5*SQU(dt/dr)/T_e0;
  std::vector<int> epsvec;
                   epsvec.resize(4); epsvec = { 1 };

  iiiprintf ( ">> efield calculation: self bias"
              " voltage used is %g, scale pot = %g,"
              " scale_efield =%g and qe=%g\n",
              Ua_SB, scale_pot, scale_efield, qe );

  for ( unsigned int r = 0; r < layer.r_dim; ++r){
    for (unsigned int z = 0; z < layer.z_dim; ++z){
      Cell& cell = layer.get_cell(r,z);
      BoundarySegment& polygon = cell.fem.boundary;

      // converting grid eps to face eps
      epsvec[0] = cell.fem.epsilon[0]; epsvec[1] = cell.fem.epsilon[1];
      epsvec[2] = cell.fem.epsilon[2]; epsvec[3] = cell.fem.epsilon[3];
      
      iiiprintf( ">> efield calculation: phi_bound un-rescaled"
                 " phi_bound[r][z]=%g, phi_bound[r+1][z]=%g,"
                 " phi_bound[r+1][z+1]=%g, phi_bound[r][z+1]=%g;"
                 " scaled with efield factor:"
                 " phi_bound[r][z]=%g, phi_bound[r+1][z]=%g,"
                 " phi_bound[r+1][z+1]=%g, phi_bound[r][z+1]=%g;"
                 " scaled with potential factor:"
                 " phi_bound[r][z]=%g, phi_bound[r+1][z]=%g,"
                 " phi_bound[r+1][z+1]=%g, phi_bound[r][z+1]=%g \n",
                 phi_bound[ r   *global_grid.mesh_z_dim+ z   ],
                 phi_bound[(r+1)*global_grid.mesh_z_dim+ z   ],
                 phi_bound[(r+1)*global_grid.mesh_z_dim+(z+1)],
                 phi_bound[ r   *global_grid.mesh_z_dim+(z+1)],
                 phi_bound[ r   *global_grid.mesh_z_dim+ z   ]/scale_efield,
                 phi_bound[(r+1)*global_grid.mesh_z_dim+ z   ]/scale_efield,
                 phi_bound[(r+1)*global_grid.mesh_z_dim+(z+1)]/scale_efield,
                 phi_bound[ r   *global_grid.mesh_z_dim+(z+1)]/scale_efield,
                 phi_bound[ r   *global_grid.mesh_z_dim+ z   ]/scale_pot,
                 phi_bound[(r+1)*global_grid.mesh_z_dim+ z   ]/scale_pot,
                 phi_bound[(r+1)*global_grid.mesh_z_dim+(z+1)]/scale_pot,
                 phi_bound[ r   *global_grid.mesh_z_dim+(z+1)]/scale_pot );

      // scaling factors for ok measures
      if (r==0) { fz = -6.0*qe; fr = -qe; } else if ( r == layer.r_dim-1 ) {
      fr = -6.0*qe/(3.0*layer.r_dim-1.0); fz = -qe/(layer.r_dim-1); } else {
      fz = -qe/r; fr = -qe/(r+1); }

      iiiprintf ( ">> efield calculation: r/z=%i/%i, "
                  "fr = %g, fz = %g scaling\n",
                  r, z, fr, fz );

      if ( std::find(epsvec.begin(),epsvec.end(),-100) != epsvec.end() ) {   // anywhere
        efield_constpot ( r, z, epsvec, Ua_SB, cell, polygon, queue ); }
      
      if ((epsvec[0] == -300 ) ||  // left
         ( epsvec[2] == -200 ) ||  // right
         ( epsvec[1] == -400 ) ) { // top
        efield_constfield ( r, z, epsvec, cell, polygon, queue ); }
      
      if ( (queue[1] == 1) && (epsvec[1] == 1) ) {
        cell.fem.efield_top = (- ( dt*cell.fem.current_top ) +
                              cell.fem.oldfield_top ); }
      if ( (queue[3] == 1) && (epsvec[3] == 1) ) {
        cell.fem.efield_bottom = (- ( dt*cell.fem.current_bottom ) +
                                 cell.fem.oldfield_bottom ); }
      if ( (queue[0] == 1) && (epsvec[0] == 1) ) {
        cell.fem.efield_left = (- ( dt*cell.fem.current_left ) +
                               cell.fem.oldfield_left ); }
      if ( (queue[2] == 1) && (epsvec[2] == 1) ) {
        cell.fem.efield_right  = (- ( dt*cell.fem.current_right ) +
                                 cell.fem.oldfield_right ); }

      if ( (z > 0) && (queue[0] == 1) && (epsvec[0] == 1) )
      { Cell& leftcell = layer.get_cell(r,z-1);
        cell.fem.efield_left += - ( dt*cell.fem.current_right );
        iiiprintf ( ">> efield calculation: leftcell found, "
                    "efield_left=%g, r=%d/%d, z=%d/%d\n",
                    cell.fem.efield_left, r, layer.r_dim-1,
                    z-1, layer.z_dim-1); } 

      if ( (r > 0) && (queue[3] == 1) && (epsvec[3] == 1) )
      { Cell& bottomcell = layer.get_cell(r-1,z);
        cell.fem.efield_bottom += - ( dt*bottomcell.fem.current_top );
        iiiprintf ( ">> efield calculation: bottomcell found, "
                    "efield_bottom=%g, r=%d/%d, z=%d/%d\n",
                    cell.fem.efield_bottom, r-1, layer.r_dim-1,
                    z, layer.z_dim-1); }

      if ( (r < layer.r_dim-1) && (queue[1] == 1) && (epsvec[1] == 1) )
      { Cell& topcell = layer.get_cell(r+1,z);
        cell.fem.efield_top += - ( dt*topcell.fem.current_bottom );
        iiiprintf(">> efield calculation: topcell found, "
                  "efield_top=%g, r=%d/%d, z=%d/%d\n",
                  cell.fem.efield_top, r+1, layer.r_dim-1,
                  z, layer.z_dim-1); }

      if ( (z < layer.z_dim-1) && (queue[2] == 1) && (epsvec[2] == 1) )
      { Cell& rightcell = layer.get_cell(r,z+1);
        cell.fem.efield_right += - ( dt*rightcell.fem.current_left );
        iiiprintf ( ">> efield calculation: rightcell found,"
                    " efield_right=%g, r=%d/%d, z=%d/%d\n",
                    cell.fem.efield_right, r, layer.r_dim-1,
                    z+1, layer.z_dim-1); }

      iiiprintf ( ">> efield calculation: r=%i/%i, z=%i/%i;"
                  " fields top=%g, bottom=%g, left=%g, right=%g; "
                  " currents top=%g, bottom=%g, left=%g, right=%g\n",
                  r, layer.r_dim-1, z, layer.z_dim-1, cell.fem.efield_top,
                  cell.fem.efield_bottom, cell.fem.efield_left,
                  cell.fem.efield_right, cell.fem.current_top,
                  cell.fem.current_bottom, cell.fem.current_left,
                  cell.fem.current_right ); 

      // reset booleans for use in loop
      queue = { 1, 1, 1, 1 };
    } // z -> layer.z_dim
  } // r -> layer.r_dim

#endif
} // void

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
