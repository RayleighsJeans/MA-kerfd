#pragma once

#include <string>
#include <vector>
#include "grid.h"
#include "boundary_check.h"
#include "output.h"
#include "var.h"

void calculate_potential( double* field_time, int r1, int r2, int z1, int z2);

void init_solver(std::string file_name, int argc, char** argv);

void init_slu(std::string file_name, std::vector<std::vector<double>> Eps); 

void destroy_solver();
void destroy_slu();

/// save phi values for const_pot and const_efield_* in phi_bound 
//
inline void boundary_conditions_to_pot(BoundaryConditionSegments bc_segments) {
  out_phi_bound("phi_bound_default.dat");
  //scale input Voltage in V to dimensionless variable 
  double scaling_factor_pot = 0.5*SQU(dt/dr)/T_e0;
  //checkout this scaling factor for constant efield 
  //(*0.5 for derivative and same scaling as for voltage)
  double scaling_factor_efield = 0.5*dr_0*scaling_factor_pot;
  for (int r = 0; r < global_grid.mesh_r_dim; ++r){
    for (int z = 0; z < global_grid.mesh_z_dim; ++z){
      for( auto&& polygon : bc_segments ){
	if(point_in_polygon(r,z,polygon) ){
	  // constant potential
	  if(polygon.type == "const_pot"  ){
	    //scale input Voltage in V to dimensionless variable 
	    phi_bound[r*global_grid.mesh_z_dim+z]    =polygon.pot*scaling_factor_pot;
	    phi_bound[(r+1)*global_grid.mesh_z_dim+z]=polygon.pot*scaling_factor_pot;
	    phi_bound[r*global_grid.mesh_z_dim+z+1]  =polygon.pot*scaling_factor_pot;
	    phi_bound[(r+1)*global_grid.mesh_z_dim+z+1]=polygon.pot*scaling_factor_pot;
	  }
	  else if(polygon.type=="const_efield_right" || polygon.type=="const_efield_left" || polygon.type=="const_efield_top"){
	    // constant electric field 
	    phi_bound[r*global_grid.mesh_z_dim+z]      =polygon.pot*scaling_factor_efield;
	    phi_bound[(r+1)*global_grid.mesh_z_dim+z]  =polygon.pot*scaling_factor_efield;
	    phi_bound[r*global_grid.mesh_z_dim+z+1]    =polygon.pot*scaling_factor_efield;
	    phi_bound[(r+1)*global_grid.mesh_z_dim+z+1]=polygon.pot*scaling_factor_efield;
	  }
	}else if(point_on_line(r,z,polygon,0.5)){
	  // constant potential
	  if(polygon.type == "const_pot"  ){
	    phi_bound[r*global_grid.mesh_z_dim+z]=polygon.pot*scaling_factor_pot;
	  }
  	  // constant electric field 
  	  else if(polygon.type=="const_efield_right" || polygon.type=="const_efield_left" || polygon.type=="const_efield_top"){
	    phi_bound[r*global_grid.mesh_z_dim+z]=polygon.efield*scaling_factor_efield;
	  }
	  // phi_bound is changed only at potential boundary, else is left at default value
	} //end: if( point_in_polygon(r,z,polygon) )
      } //end: for( auto&& polygon : bc_segments )
    }
  }
  out_phi_bound("phi_bound.dat");
}


/// fill espilon and phi_bound
inline void set_solver_boundary(std::string file_name, std::vector<std::vector<double>>& Eps) 
{
  //auto bc_segments = read_boundary_conditions(file_name);
  
  // the new epsilon filling function
  init_bounds(file_name, Eps);
  // printing eps matrix on rank 0
  out_eps(Eps);
  // printing phi boundaries on rank 0
  out_phi_bound("phi_bound.dat");

}
