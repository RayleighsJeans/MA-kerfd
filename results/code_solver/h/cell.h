#pragma once

#include "pic.h"
#include "var_diag.h"
#include <stdlib.h>
#include <vector>
#include <array>
#include <iostream>

// TODO this might be better in the wiki
// This approach is used to keep the particles 
// array linear and without gaps.
// If it would contain gaps these would need to be closed
// after the pusher is done.
// Keep in mind that order inside the particle array does not matter.
// 
// +-+-+-+-+-+-+-+-+-+-+-+-+
// |x|x|x|x|x|x|x|x|x|x|x|x|
// +-+-+-+-+-+-+-+-+-+-+-+-+
//  ^                       ^
//  0                       1
//                          2
//                          3
// 
// 0 = begin of the particle array
// 1 = end particles of the particle array
// 2 = begin of the new particles that will arrive
// 3 = end of the new particles that will arrive
// 
// Case 1.0: new particles arrive 
// 
// +-+-+-+-+-+-+-+-+-+-+-+-+-+
// |x|x|x|x|x|x|x|x|x|x|x|x|y|
// +-+-+-+-+-+-+-+-+-+-+-+-+-+
//  ^                       ^ ^
//  0                       1
//                          2
//                            3
// 
// new particles are now stored from 2 to exclusive 3
// old particle between 0 and exclusive 1
// 
// Case 1.1: new particles arrive and old are deleted
// 
// +-+-+-+-+-+-+-+-+-+-+-+-+-+
// |x|x|x|x|x|x|x|x|x|x|x|x|y|
// +-+-+-+-+-+-+-+-+-+-+-+-+-+
//  ^         ^             ^ ^
//  0         5             1
//                          2
//                            3
// 5 = old particle to be deleted
// 
// if the particle that will be used to fill in the gap of 5 
// has an lower than 2 this means that you fill the gap with a
// particle that has to be evaluated.
// if it is no lower than 2 it means that it is an incoming particle
// and has to be skipped because it was already pushed by another cell
// to this cell
// 
// +-+-+-+-+-+-+-+-+-+-+-+-+
// |x|x|x|x|x|y|x|x|x|x|x|x|
// +-+-+-+-+-+-+-+-+-+-+-+-+
//  ^         ^             ^
//  0         6             1
//                          2
//                          3
// 
// 6 = new particle is stored inside the normal array

//Boundary Types (the ones with direction are used for domain boundary conditions, ORDER OF DIRECTION IS CRUCIAL! (LEFT->TOP->RIGHT))
enum BoundaryType{
  NO_BOUNDARY,
  OUT_OF_DOMAIN,
  DIEL_BORDER,
  DIELECTRIC,
  REFLECTION,
  REEMISSION,
  RECYCLE,
  SEE_DIEL,
  MAX_BOUNDARY_TYPES,	//8

  LEFT_NO_BOUNDARY,
  LEFT_OUT_OF_DOMAIN,
  LEFT_DIEL_BORDER,
  LEFT_DIELECTRIC,
  LEFT_REFLECTION,
  LEFT_REEMISSION,
  LEFT_RECYCLE,
  LEFT_SEE_DIEL,
  MAX_BOUNDARY_TYPES_LEFT,  //17

  TOP_NO_BOUNDARY,
  TOP_OUT_OF_DOMAIN,
  TOP_DIEL_BORDER,
  TOP_DIELECTRIC,
  TOP_REFLECTION,
  TOP_REEMISSION,
  TOP_RECYCLE,
  TOP_SEE_DIEL,
  MAX_BOUNDARY_TYPES_TOP,  //26

  RIGHT_NO_BOUNDARY,
  RIGHT_OUT_OF_DOMAIN,
  RIGHT_DIEL_BORDER,
  RIGHT_DIELECTRIC,
  RIGHT_REFLECTION,
  RIGHT_REEMISSION,
  RIGHT_RECYCLE,
  RIGHT_SEE_DIEL,
  MAX_BOUNDARY_TYPES_RIGHT //35
};

struct Diagnostic{
 unsigned int particle_number;     //particle number per cell
#if DIAG_FLUX
 std:: vector<int> particle_flux;  //particle flux per cell, positive in x and z direction
 std:: vector<double> energy_flux;    //energy flux per cell
#endif

#if DIAG_COLL
 unsigned int coll_ioniz;	    //counts ionization collisions
#endif
#if DIAG_VELZ
 std:: vector<double> particle_velocity; 
 std:: vector<double> particle_velocity_squared;
#endif

#if DIAG_DIST
 std:: vector<int> vel_r_dist;
 std:: vector<int> vel_z_dist;
 std:: vector<int> vel_t_dist; //TODO check if this is really used
 std:: vector<int> vel_dist;
#endif

 //Phase resolved diagnostics for RF
#if DIAG_PR
 std::vector<int> particle_number_pr;
 std::vector<int> pr_count;
#if DIAG_VELZ
 std::vector<double> velzr_pr;
 std::vector<double> velzz_pr;
 std::vector<double> velzt_pr;
#endif
#if DIAG_DIST
 std::vector<std::vector<int>> vel_r_dist_pr;
 std::vector<std::vector<int>> vel_z_dist_pr;
 std::vector<std::vector<int>> vel_t_dist_pr;
 std::vector<std::vector<int>> vel_dist_pr;
#endif
#if DIAG_COLL
 std::vector<int> coll_ioniz_pr;
#endif
#endif

 Diagnostic(){
#if DIAG_VELZ
    particle_velocity.resize(3);
    particle_velocity_squared.resize(3);
#endif

    //particle flux (up, right, down, left)
#if DIAG_FLUX
    particle_flux.resize(4);
    energy_flux.resize(4); 
#endif

#if DIAG_DIST
    //TODO make this flexible to cells with no particle in them
    //vel_r_dist.resize(NBIN);  
    //vel_z_dist.resize(NBIN);
    //vel_t_dist.resize(NBIN);
    //vel_dist.resize(NBIN);
#if DIAG_PR
    //TODO make this flexible as well, otherwise you will run out of memory
    vel_r_dist_pr.resize(NBIN);  
    vel_z_dist_pr.resize(NBIN);
    vel_dist_pr.resize(NBIN);
#endif
#endif
}

 void clear(){
      particle_number = 0;
#if DIAG_VELZ
      for(int i = 0; i < particle_velocity.size(); i++){
	particle_velocity[i]=0.0;
	particle_velocity_squared[i]=0.0;
      }	
#endif
#if DIAG_FLUX
      for(int i = 0; i < particle_flux.size(); i++){
	//particle_flux[i]=0;	
	energy_flux[i]=0.0;
      }	
#endif
#if DIAG_COLL
      coll_ioniz =0;
#endif
#if DIAG_DIST
      for(int i=0; i<vel_r_dist.size(); i++){
	vel_r_dist[i]=0;
	vel_z_dist[i]=0;
	vel_t_dist[i]=0;
	vel_dist[i]=0;
      }
#endif
 }

};

struct fem_diag{
#if USE_FEM_SOLVER
  // area weighted charge for fem
  double area_weighted_charge;
  // calculated current for fem
  double current_top;
  double current_bottom;
  double current_right;
  double current_left;
  // efield from fem
  double efield_top;
  double efield_bottom;
  double efield_right;
  double efield_left;
  // oldfield from fem
  double oldfield_top;
  double oldfield_bottom;
  double oldfield_right;
  double oldfield_left;
#endif

  void clear_fem (){
  #if USE_FEM_SOLVER
    // area charge
    area_weighted_charge = 0.0;
    // current
    current_top = 0.0;
    current_bottom = 0.0;
    current_left = 0.0;
    current_right = 0.0;
    // efield
    efield_top = 0.0;
    efield_bottom = 0.0;
    efield_left = 0.0;
    efield_right = 0.0;
  #endif
  }
  
  void clear_oldfield(){
  #if USE_FEM_SOLVER
    // oldfield
    oldfield_top = 0.0;
    oldfield_bottom = 0.0;
    oldfield_left = 0.0;
    oldfield_right = 0.0;
  #endif
  }

};

//Struct which contains boundary properties of the eight Neighbourcells
struct NeighbourCellProps {
  //direction of bounds is starting from the upper neighbour and then counter clockwise
  std::vector<int> bounds;  
  //if bound is either REFLECTON or REEMISSION ior SEE save normal vectors 
  std::array<std::array<double,2>,8> normals = {{{{0,0}},{{0,0}},{{0,0}},{{0,0}},{{0,0}},{{0,0}},{{0,0}},{{0,0}}}};
  // 1---0---7
  // |	     |
  // 2       6
  // |       |
  // 3---4---5
};

struct Cell {

  Cell( ) {
    push_mode=false;
    dens.fill(0);
    surf_dens.fill(0);
    recycle.fill(0);
    boundary_type=NO_BOUNDARY;
  }

  void push_back( Particle& pt ) {
    particles.push_back( pt );
  }

  void emplace_back( Particle&& pt ) {
    particles.emplace_back( pt );
  }

  Particle remove( unsigned int& id ){
    Particle ret = particles[id];
    particles[id] = particles.back(); 
    particles.pop_back();
    if ( (particles.size()  < end_particles) && (push_mode==true) ){
      id--;
      end_particles--;
    }else if(push_mode==false){
      id--;
      end_particles=particles.size();
    }
    return ret;
  }

  size_t size(){
    return particles.size();
  }

  bool empty(){
    return particles.empty();
  }

  void clear(){
    particles.clear();
  }
  
  void shrink_to_fit(){
    particles.shrink_to_fit();
  }

  size_t memory_consumption(){
    return particles.capacity() * sizeof(Particle);
  }

  void register_prepush_state(){
    push_mode=true;
    end_particles = particles.size();
  }

  void register_postpush_state(){
    push_mode=false;
    end_particles = particles.size();
  }

  int sanity_check( unsigned int r, unsigned int z ){
    int dont_belong_here = 0;
    for( auto& particle : particles ){
      unsigned int r_p = particle.r;
      unsigned int z_p = particle.z;
      if ( (r_p != r) || (z_p != z) ) {
	dont_belong_here++;
	//std::cout<<" particle at "<<particle.r<<","<<particle.z<<" is outside of cell "<<r<<","<<z<<std::endl;
      }
    }
    return dont_belong_here;
  }

  void print(){
    for (int i = 0; i < particles.size(); ++i){
      std::cout << particles[i].r << " ";
    }
    std::cout << std::endl;
  }

  // the vector holding all particle information
  std::vector<Particle> particles; 

  // needed for book keeping
  size_t end_particles;
  bool push_mode;


  fem_diag fem;
  
  void clear_oldfield(){
    fem.clear_oldfield();
  }
  
  void clear_fem(){
    fem.clear_fem();
  }

  Diagnostic diagnostic_arrays;
  
  void clear_diagnostics(){
    diagnostic_arrays.clear();
  }

  // interpolated particle density for solver (0= up left and counter clockwise)
  std::array<double,4> dens;
  // 0-------3
  // |	     |
  // |       |
  // |       |
  // 1-------2
  //surface density vector for all four corners of cell (0=up left, 1=down left, 2=down right, 3= up right)+neighbours
  std::array<double,12> surf_dens;
  //         4	     11
  //         |	     |
  //         |	     |
  //         |	     |
  // 5-------0-------3-------10
  //	     |	     |
  //	     |       |
  //	     |       |
  // 6-------1-------2-------9
  //         |	     |
  //         |	     |
  //         |	     |
  //         7	     8
  
  //counts the particles impinging on surfaces for future reinjection of particles of a different species(0=up, 1= left, 2=down, 3=right)
  std::array<int,4> recycle;
  //  ---0---
  // |	     |
  // 1       3
  // |       |
  //  ---2---
  
  // properties needed for boundary conditions
  int r = -1;
  int z = -1;
  
  bool touched_by_polygon = false;
  bool boundary_cell = false;
  enum BoundaryType boundary_type; 

  //normalized normal vector of boundary condition given, for NO_BOUNDARY simply (0,0), for the other cells (r_norm,z_norm)
  std::vector<double> normal_vector;
  
  bool is_dielectric;
  bool is_dielectric_surface;
  
  NeighbourCellProps neighbours;

};

