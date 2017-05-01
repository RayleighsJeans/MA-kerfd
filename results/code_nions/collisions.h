#pragma once

// for particles 
#include "pic.h"
// for reactions
#include "react.h"
// for grid and gridlayer
#include "grid.h"

 void coll_coulomb( GridLayer& layer, double Acoll, Vec3d *momcheck, double *engcheck, 
	 int ncoll_coulomb[], double coll_dr[], double coll_dt[] );
 

   
 void coll_el_all_fake(    GridLayer& neutral_layer,		   // molecules
                           GridLayer& electron_layer,              // electrons
                           Reaction React,
                           Vec3d *momcheck, double *engcheck, 
			   int ncoll_el_ntrl[], double ecoll_dr[], double ecoll_dt[], double ncoll_dr[], double ncoll_dt[] );
                           

                          
 void coll_el_ntrl_ioniz( GridLayer& neutral_layer, 
                          GridLayer& electron_layer,
                          GridLayer& ion_layer,
			  Reaction React,
                          Vec3d *momcheck, double *engcheck, int ncoll_ioniz[], 
			  double ecoll_dr[], double ecoll_dt[], double ncoll_dr[], double ncoll_dt[] );
                                                    
                                
                          
void coll_ion_ntrl(      GridLayer& neutral_layer,		  // molecules
                         GridLayer& ion_layer,			  // ions
                         Reaction React, int ncoll_ion_ntrl[],
			 double icoll_dr[], double icoll_dt[], double ncoll_dr[], double ncoll_dt[] );  
                        
void coll_ntrl_ntrl_acc(   Particle  neutrals[],   int ordcount_n[],  int nr, int nz,     Reaction React,  Vec3d *momcheck, 
			double *engcheck, int ncoll_ntrl_ntrl[], double ncoll_dr[], double ncoll_dt[] );

void coll_ntrl_ntrl_weighted(   Particle  neutrals[],   int ordcount_n[],  int nr, int nz,     Reaction React,  Vec3d *momcheck, 
			double *engcheck, int ncoll_ntrl_ntrl[], double ncoll_dr[], double ncoll_dt[] );


void coll_ntrl_ntrl( GridLayer& layer,
		    Reaction React, Vec3d *momcheck, double *engcheck,
		    int ncoll_ntrl_ntrl[], double ncoll_dr[], double ncoll_dt[]
    );

 void coll_el_ntrl_ioniz2( GridLayer& neutral_layer, GridLayer& electron_layer, GridLayer& ion_layer,
			  Reaction React, Vec3d *momcheck, double *engcheck, int ncoll_double_ioniz[],
                        double ecoll_dr[],double ecoll_dt[], double ncoll_dr[],double ncoll_dt[]);


void calculate_pairs_per_cell( int ordcount_n[], Reaction react, int nr, int nz, long pairs_per_cell[]);

void coll_el_diss_attach( GridLayer& neutral_layer, // was dt_ion,     // Check the ion velocity!!
          GridLayer& electron_layer,
          GridLayer& nion_layer, 
	//  GridLayer& atom_layer
	  double ma_over_me,
          Reaction React_O2_dat, Vec3d *momcheck, double *engcheck, 
	  int ncoll_ioniz[], double ecoll_dr[], double ecoll_dt[], double ncoll_dr[], double ncoll_dt[]  ); 

void coll_el_diss_recomb_fake(GridLayer& neutral_layer, double M_n, int dt_ntrl,
                          GridLayer& electron_layer,
			  Reaction React,
                          Vec3d *momcheck, double *engcheck, int ncoll_ioniz[],
			  double ecoll_dr[], double ecoll_dt[], double ncoll_dr[], double ncoll_dt[] );

 void coll_el_detachment( Particle  neutrals[],  int *nn,  int ordcount_ntrls[], double M_n, int dt_ntrl,
                          Particle  electrons[], int *ne,  int ordcount_el[],
                          Particle  ions[],      int *ni,  
			  int nr,   int nz,  Reaction React,
                          Vec3d *momcheck, double *engcheck, int ncoll_ioniz[],
			  double ecoll_dr[], double ecoll_dt[], double ncoll_dr[], double ncoll_dt[] );

 void coll_diss_ionization(  Particle  molecules[],  int *nm,  int ordcount_m[], double M_m,   // molecules
                             Particle  electrons[],  int *ne,  int ordcount_el[],              // electrons
                             Particle  ions[],       int *ni,  double M_i,
                             Particle  neutrals[],   int *nn,  double M_n,  int H_out,
                             int nx,    Reaction React,
                             Vec3d *momcheck, double *engcheck   );




 void coll_ip_im_neutralization_fake( GridLayer&  neutral_layer,
                          GridLayer& ion_layer, 
			  Reaction React,
                          Vec3d *momcheck, double *engcheck, int ncoll_ioniz[],
			  double ecoll_dr[], double ecoll_dt[], double ncoll_dr[], double ncoll_dt[] );

 void coll_neutr_detachment_fake(GridLayer& neutral_layer,
                          GridLayer& ion_layer,
                          GridLayer& electron_layer,
			  Reaction React,
                          Vec3d *momcheck, double *engcheck, int ncoll_ioniz[],
			  double ecoll_dr[], double ecoll_dt[], double ncoll_dr[], double ncoll_dt[] );



void collisions( double* engcheck,
		    Vec3d* momcheck,
		    double* coll_time,
		    double* order_time);

void neutral_collisions(
		    double* engcheck,
		    Vec3d* momcheck,
		    double* colN_time,
		    double* orderN_time,
		    Reaction reaction
);
