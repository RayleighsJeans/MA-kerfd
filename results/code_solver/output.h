#pragma once

#include <vector>
#include <time.h>
#include <string>
#include "var_diag.h"
#include "init_diag.h"
#include "grid_layer.h"
#include "boundary_check.h"
#include "grid.h"

void out_test();

void out_dens( GridLayer& layer, std::string dat_name );

void out_dens_pr( GridLayer& layer, std::string dat_name );

void out_and_zero_phi( double phi_aver[], double dt, double dr, std::string dat1 );

void out_phi_bound();

void out_and_zero_coll( int ncoll[], std::string dat1 );

void out_flux_P( GridLayer& layer, int nr_model, int z_pos, std::string dat_i);

void out_flux_E( GridLayer& layer, int nr_model, int z_pos, std::string dat_i);

void out_flux_P_sum( GridLayer& el_layer, GridLayer& ion_layer, GridLayer& ntrl_layer, int nr_model, int z_pos, std::string dat_sum);

void out_flux_E_sum( GridLayer& el_layer, GridLayer& ion_layer, GridLayer& ntrl_layer, int nr_model, int z_pos, std::string dat_sum);

void out_velz( GridLayer& layer, std::string dat1, std::string dat2, std::string dat3 );

void out_velz_pr(GridLayer& layer, std::string datr, std::string datz, std::string datt);

void out_tempz( GridLayer& layer, std::string dat1, std::string dat2 );

void out_debye(GridLayer &layer, std::string dat1, std::string dat2, std::string dat3);

void out_vel_fluidlines(GridLayer& layer, int n, std::string dat);

void out_flux(GridLayer& layer, std::string dat1, std::string dat2, std::string dat3, std::string dat4);

void out_dist(GridLayer& layer, std::string fnDistz, std::string fnDistr, std::string fnDist);

void out_and_zero_dist_pr(GridLayer& layer, std::string fnDistz, std::string fnDistr, std::string fnDist);

void out_mfp(double ecoll_dr[],double ecoll_dt[],double icoll_dr[],double icoll_dt[],double ncoll_dr[],double ncoll_dt[], std::string Fcoll,int nBin);

void out_coll_pr(GridLayer& layer, std::string dat);

void out_ang_curr_enDist( long int& ions_tot, std::vector<double>& Angle, std::vector<int>& Bin_particles, 
			  std::vector<double>& Energy, std::vector<std::vector<int>>& Bin_particle_energy, 
			  std::vector<std::vector<std::vector<int>>>& Bin_particle_origin,
			  double nav_steps, std::string  f_name );

void out_eps(std::vector<std::vector<double> >Eps);
void out_bounds(GridLayer& layer, std::string f);
void out_recycle(GridLayer& layer, std::string f);
void out_neighbour_bounds(GridLayer& layer, std::string f);
void out_normals(GridLayer& layer, std::string f);
void out_neighbour_normals(GridLayer& layer, std::string f);
void out_phi_bound(std::string f);
void out_rhs( double phi_av[], double dt, double dr, std::string dat1 );
void out_efield_dc(std::vector<double> efield_dc);
void out_current_self_bias(std::vector<double> efield_dc, std::vector<double> I_el, std::vector<double> I_ion);
void out_phi_scaled(double* phi);
void out_efield_scaled(Field* efield);

void out_boundary_conditions_to_gnuplot(Grid& grid, BoundaryConditionSegments& bc_segments);
void out_variables_to_gnuplot();

void out_and_zero_current_to_walls( std::vector<CountParticles>& count_particles);

void out_surf_dens(GridLayer& layer, std::string file_name);

void out_pt_too_fast(GridLayer& layer, std::string file_name);

void out_charge_dens(GridLayer& layer, std::string file_name);

  void write_diagnostics( 
    double field_time, 
    double pushE_time, 
    double coll_time, 
    double initN_time, 
    double pushN_time,
    double colN_time,
    double inj_time,
    double engcheck, 
    Bench_Start begin_main_time,
    Vec3d momcheck,
    std::string datatrc
);

void write_diagnostics_cell_based();

void reduce_for_output();
