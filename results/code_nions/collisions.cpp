
#include "collisions.h"
#include "init_prtcls.h"
#include "debug_printing.h"
#include "domain_decomposition.h"
	  
static inline void get_distinct_particle_ids( int& n1, int& n2, int max ){
  assert( max >= 2 && " max must be at least 2 in order to return 2 distinct elements " );
  // this correctos this to never generate the max value making it a half open range 
  // between 0 and max == [ 0, max )
  auto corrector = 1.0 - std::numeric_limits<double>::min();
  n1 = RAND * corrector * max;
  n2 = RAND * corrector * max;
  if ( n1 == n2 ) {
    // cannot go up
    if ( n2 + 1 == max ){
      n2--;
    } 
    // cannot go down
    if ( n2 - 1 == -1 ){
      n2++;
    } 
  }
}


static inline void diag_collisions( Particle& p, double* coll_dr, double* coll_dt ) 
{
#if DIAG_COLL
  // mean free path of a species
  if ((int)(p.coll_r / d_mfp) >= 0 &&
      (int)(p.coll_r / d_mfp) < n_bins_coll) {
    coll_dr[(int)(p.coll_r / d_mfp)] += 1.;
  } else {
    coll_dr[ n_bins_coll - 1] += 1.;
  }

  // mean free time of a species //do not calculate if restart from a run with no coll diag, beacuse p.coll_t will be zero
  // simply avoid the diagnostic till particle collided once
  if(p.coll_t!=0){
    if( (int)(nstep - p.coll_t) >= 0 &&
	(int)( (nstep - p.coll_t) / d_mft) < n_bins_coll ) {
      //printf("nstep - p.coll_t = %e -> bin %i\n",nstep - p.coll_t, (int)( (nstep - p.coll_t) / d_mft));
      coll_dt[(int)((nstep - p.coll_t) / d_mft)] += 1.;
    } else {
      coll_dt[ n_bins_coll - 1] += 1.;
    }
  }
  p.coll_r = 0.;
  p.coll_t = nstep;
#endif
}


/***********************************************************************************************************
e-e or i-i Coulomb collision routine (for the same type of particles)
kind = 0 e-e
otherwise i-i
xyz -> rtz
************************************************************************************************************/

void coll_coulomb( GridLayer& layer,
		  double Acoll, Vec3d *momcheck, double *engcheck,
		  int ncoll_coulomb[], double coll_dr[], double coll_dt[]) {
  Vec3d v_rel;
  double W2, W, IW, v_proj_rt, ivp, cos_theta, sin_theta, cos_beta, sin_beta;
  double cos_phi, sin_phi, psi;
  double vr, vt, vz;
  double tg2_phim_d2;
  int N2coll, loopind, pind;
  int inds2coll[NEcell];

  double VolFact;

  for ( int i = 0; i < layer.r_dim; ++i) {
    // Zylindrical volume scaling according  Verboncoer,JCP 2001, equ: 19 + 20
    if (i != 0)
      VolFact = 1. / i;
    else
      VolFact = 6.;

    for ( int j = 0; j < layer.z_dim; ++j) {
      auto& cell = layer.get_cell( i , j );
      N2coll = cell.size();

      if (N2coll > NEcell)
	fprintf(stderr,
		"To much electrons !!! increase NEcell ir=%d jz=%d  n_e=%zu", i,
		j, cell.size());

      if (N2coll > 1) {
	for (loopind = 0; loopind < N2coll; ++loopind)
	  inds2coll[loopind] = loopind;  // nachalniy massiv indexov chastic

	int k, m;
	while (N2coll > 0) {
	  // TODO check wheather m has to be calculated in the else branch
	  if (N2coll > 1) {
	    pind = (int)(N2coll * RAND);
	    if (pind == N2coll) --pind;  // a mozhet ne nuzhno ?
	    m = inds2coll[pind];
	    for (loopind = pind + 1; loopind < N2coll; ++loopind)
	      inds2coll[loopind - 1] = inds2coll[loopind];
	    --N2coll;

	    pind = (int)(N2coll * RAND);
	    if (pind == N2coll) --pind;  // a mozhet ne nuzhno ?
	    k = inds2coll[pind];
	    for (loopind = pind + 1; loopind < N2coll; ++loopind)
	      inds2coll[loopind - 1] = inds2coll[loopind];
	    --N2coll;
	  } else {
	    k = inds2coll[N2coll - 1];
	    --N2coll;
	  }

	  // Just for diagnostic purposes - remove it when U done BEGIN
	  
	  Particle& p1 = cell.particles[m];
	  Particle& p2 = cell.particles[k];

	  (*momcheck).r += p1.vr + p2.vr;
	  (*momcheck).t += p1.vt + p2.vt;
	  (*momcheck).z += p1.vz + p2.vz;

	  *engcheck += SQU(p1.vr) + SQU(p1.vt) + SQU(p1.vz) 
	             + SQU(p2.vr) + SQU(p2.vt) + SQU(p2.vz);

	  v_rel.r = p1.vr - p2.vr;  // relative velocity
	  v_rel.t = p1.vt - p2.vt;
	  v_rel.z = p1.vz - p2.vz;

	  W2 = SQU(v_rel.r) + SQU(v_rel.t) + SQU(v_rel.z);
	  W = sqrt(W2);	// Absolutnaya velichina otn. skorosti
	  W = MAX(W, 1.e-10);  // A chtob na 0 ne delit'
	  IW = 1. / W;

	  v_proj_rt = sqrt(SQU(v_rel.r) + SQU(v_rel.t));  // proekciya na XY
	  v_proj_rt = MAX(v_proj_rt, 1.e-10);		  // ta zhe fignia
	  ivp = 1. / v_proj_rt;

	  cos_theta = v_rel.z * IW;  // theta - ugol mezhdu v_rel i OZ
	  sin_theta = v_proj_rt * IW;

	  cos_beta =
	      v_rel.r *
	      ivp;  // beta - ugol mezhdu proekciey v_rel na XY i os'yu abciss
	  sin_beta = v_rel.t * ivp;

	  tg2_phim_d2 = -2. * log(RAND + 1.e-20);
	  cos_phi = cos(TWOPI * RAND);  // not good idea

	  tg2_phim_d2 *= SQU(cos_phi) * Acoll * cell.size() * IW * IW * IW *
			 VolFact;  // last vesrsion  - Gauss

	  cos_phi = (1. - tg2_phim_d2) / (1. + tg2_phim_d2);

	  if (cos_phi < -1. || cos_phi > 1.) cos_phi = 1. - 2. * RAND;

	  sin_phi = sqrt(1. - SQU(cos_phi));
	  psi = TWOPI * RAND;

	  vz = W * cos_phi;  //   komponenty posle stolknoveniya v povernutoy
			     //   coor. sys.
	  vr = W * sin_phi * cos(psi);
	  vt = W * sin_phi * sin(psi);

	  v_rel.r -= vr * cos_beta * cos_theta - vt * sin_beta +
		     vz * cos_beta *
			 sin_theta;  // eto  uzhe minus prirascheniya otn. skor.
	  v_rel.t -= vr * sin_beta * cos_theta + vt * cos_beta +
		     vz * sin_beta * sin_theta;  // a ne otn. skorost'
	  v_rel.z -= vz * cos_theta - vr * sin_theta;

	  p1.vr -= 0.5 * v_rel.r;  // Note inversed sign due inversed delta be4
	  p1.vt -= 0.5 * v_rel.t;
	  p1.vz -= 0.5 * v_rel.z;

	  p2.vr += 0.5 * v_rel.r;
	  p2.vt += 0.5 * v_rel.t;
	  p2.vz += 0.5 * v_rel.z;

// collision diagnostics
#if DIAG_COLL
	 ncoll_coulomb[i*layer.z_dim+j]++;
	 diag_collisions( p1, coll_dr, coll_dt );
	 diag_collisions( p2, coll_dr, coll_dt );
#endif

	  // Just for diagnostic purposes - remove it when U done **BEGIN

	  (*momcheck).r -= p1.vr + p2.vr;
	  (*momcheck).t -= p1.vt + p2.vt;
	  (*momcheck).z -= p1.vz + p2.vz;

	  *engcheck -= SQU(p1.vr) + SQU(p2.vr) +
		       SQU(p1.vz) + SQU(p2.vz) +
		       SQU(p1.vt) + SQU(p2.vt);

	  //                                                      **END (of
	  //                                                      diagnostics )
	}
      }
    }
  }
}

/************************************************************************************************
  Simplified version of inelastic collisions of electron with whatsoever -
  we just let electron loose some energy and then rotate it's velocity
  we don't give a damn to reaction products!
*************************************************************************************************/
void coll_el_all_fake(GridLayer& neutral_layer, // molecules
		      GridLayer& electron_layer,  // electrons
		      Reaction React, Vec3d *momcheck,
		      double *engcheck, int ncoll_el_ntrl[], double ecoll_dr[],
		      double ecoll_dt[], double ncoll_dr[], double ncoll_dt[]) {
  const double m_e = 1.;

  double dtn, VolFact;
  double W2, W, IW, W2_0, W_0, IW_0, v_proj_rt, ivp, cos_theta, sin_theta,
      cos_beta, sin_beta;
  double cos_phi, sin_phi, psi;

  double vr, vt, vz;
  Vec3d v_rel, delta_v_rel, target_v;

  const double M_m=neutral_layer.mass;
  const int dt_ntrl=neutral_layer.dt_subcycling;

  int count_coll = 0;

  double E_th, Emin, Emax, Estep;
  int Epoints, Eind;

  E_th = React.Eth;
  Emin = React.Emin;
  Emax = React.Emax;
  Estep = React.Estep;
  Epoints = React.N;

  dtn = 1. / dt_ntrl;  // Scaling factor 4 ntrl velocity 

  for (int i = 0; i < neutral_layer.r_dim; ++i) {
    // Zylindrical volume scaling according  Verboncoer,JCP 2001, equ: 19 + 20
    if (i != 0)
      VolFact = 1. / i;
    else
      VolFact = 6.;

    for (int j = 0; j < neutral_layer.z_dim; ++j) {
      auto& neutral_cell = neutral_layer.get_cell( i, j );
      auto& electron_cell = electron_layer.get_cell( i, j );

      if ( !neutral_cell.empty() && !electron_cell.empty() ) {

	for (int i_el = 0; i_el < electron_cell.size(); i_el++) {
	  int i_m = (int)(neutral_cell.size() * RAND);
	  if (i_m == neutral_cell.size()) i_m--;

	  Particle& neutral = neutral_cell.particles[i_m];
	  Particle& electron = electron_cell.particles[i_el];

	  count_coll++;
	  v_rel.r = electron.vr - neutral.vr * dtn;  // relative velocity
	  v_rel.t = electron.vt - neutral.vt * dtn;  // !!3.09.01 Check it on scaling 4 ion velocity
	  v_rel.z = electron.vz - neutral.vz * dtn;

	  W2_0 =
	      SQU(v_rel.r) + SQU(v_rel.t) + SQU(v_rel.z);  // before scattering
	  W_0 = sqrt(W2_0);

	  // Linear fit for Cross-section
	  double S_i = 0.;
	  if (W2_0 >= E_th) {
	    Eind = (int)((W2_0 - Emin) / Estep);
	    if (Eind < 0)
	      S_i = React.CS[0];
	    else if (Eind >= Epoints)
	      S_i = React.CS[Epoints];
	    else
	      S_i = React.CS[Eind] + (React.CS[Eind + 1] - React.CS[Eind])*((W2_0 - Emin) / Estep - Eind);
	  }

	  
	  if ( W2_0 >= E_th && RAND < neutral_cell.size() *W_0*S_i*VolFact ) //!! Est' collision!! 1-exp(-n_m*W_0*Si*dt_coll)
	  {
	    // Start of real collision routine

	    *engcheck += 0.5 * m_e * (SQU(electron.vr) +
				      SQU(electron.vt) +
				      SQU(electron.vz));
	    *engcheck += 0.5 * M_m * (SQU(neutral.vr * dtn) +
				      SQU(neutral.vt * dtn) +
				      SQU(neutral.vz * dtn)) -
			 0.5 * m_e * M_m / (m_e + M_m) * E_th;

	    (*momcheck).r += m_e * electron.vr + M_m * neutral.vr * dtn;
	    (*momcheck).t += m_e * electron.vt + M_m * neutral.vt * dtn;
	    (*momcheck).z += m_e * electron.vz + M_m * neutral.vz * dtn;

	    W_0 = MAX(W_0, 1.e-14);
	    IW_0 = 1. / W_0;

	    W2 = W2_0 - E_th;
	    W = sqrt(W2);  // Absolutnaya velichina otn. skorosti

	    delta_v_rel.r = v_rel.r * (W * IW_0 - 1.);  // After energy loss
	    delta_v_rel.t = v_rel.t * (W * IW_0 - 1.);
	    delta_v_rel.z = v_rel.z * (W * IW_0 - 1.);

	    v_rel.r += delta_v_rel.r;
	    v_rel.t += delta_v_rel.t;
	    v_rel.z += delta_v_rel.z;

	    W = MAX(W, 1.e-10);  // A chtob na 0 ne delit'
	    IW = 1. / W;

	    v_proj_rt = sqrt(SQU(v_rel.r) + SQU(v_rel.t));  // proekciya na XY
	    v_proj_rt = MAX(v_proj_rt, 1.e-10);		    // ta zhe fignia
	    ivp = 1. / v_proj_rt;

	    cos_theta = v_rel.z * IW;  // theta - ugol mezhdu v_rel i OZ
	    sin_theta = v_proj_rt * IW;

	    cos_beta = v_rel.r * ivp;  // beta - ugol mezhdu proekciey v_rel na XY i os'yu abciss
	    sin_beta = v_rel.t * ivp;

	    /* Sampling scattering angles in C.M. system */
	    cos_phi = 2. * RAND - 1.;  //FIXME: just 4 test!! Use = (2.+E-2.*pow((1+E), RAND))/E instead
	    sin_phi = sqrt(1. - SQU(cos_phi));
	    psi = TWOPI * RAND;

	    vz = W * cos_phi;  //   komponenty posle stolknoveniya v povernutoy
			       //   coor. sys.
	    vr = W * sin_phi * cos(psi);
	    vt = W * sin_phi * sin(psi);

	    v_rel.r -= vr*cos_beta*cos_theta - vt*sin_beta + vz*cos_beta*sin_theta;  // eto  uzhe minus prirascheniya otn. skor.
	    v_rel.t -= vr*sin_beta*cos_theta + vt*cos_beta + vz*sin_beta*sin_theta;  // a ne otn. skorost'
	    v_rel.z -= vz*cos_theta - vr*sin_theta;

	    v_rel.r -= delta_v_rel.r;
	    v_rel.t -= delta_v_rel.t;
	    v_rel.z -= delta_v_rel.z;

	    target_v.r = neutral.vr * dtn;  // remove it after debuging
	    target_v.t = neutral.vt * dtn;  // remove it after debuging
	    target_v.z = neutral.vz * dtn;  // remove it after debuging

	    electron.vr -= M_m / (M_m + m_e) * v_rel.r;
	    electron.vt -= M_m / (M_m + m_e) * v_rel.t;
	    electron.vz -= M_m / (M_m + m_e) * v_rel.z;

	    target_v.r += m_e / (M_m + m_e) * v_rel.r;
	    target_v.t += m_e / (M_m + m_e) * v_rel.t;
	    target_v.z += m_e / (M_m + m_e) * v_rel.z;

	    *engcheck -= 0.5 * m_e * (SQU(electron.vr) +
				      SQU(electron.vt) +
				      SQU(electron.vz));
	    *engcheck -= 0.5 * M_m * (SQU(target_v.r) + SQU(target_v.t) + SQU(target_v.z));

	    (*momcheck).r -= m_e * electron.vr + M_m * target_v.r;
	    (*momcheck).t -= m_e * electron.vt + M_m * target_v.t;
	    (*momcheck).z -= m_e * electron.vz + M_m * target_v.z;

// collision diagnostics
#if DIAG_COLL
	    ncoll_el_ntrl[i*neutral_layer.z_dim+j]++;
	    diag_collisions( electron, ecoll_dr, ecoll_dt );
	    diag_collisions( neutral, ncoll_dr, ncoll_dt );
#endif

	  }  //  if (W2_0 >= E_th && RAND < n_m*W_0*Si)
	}  // for ( i_el = 0; i_el < n_el; i_el++)
      }  // if ( n_m && n_el )
    }  // for (i=0; i < nr;  i++)
  }
}

/************************************************************************************************
  Simplified version of elastic collisions of ion with whatsoever -
  we just let ion loose some energy and then rotate it's welocity
  we don't give a damn to reaction products!
*************************************************************************************************/
void coll_ion_ntrl(GridLayer& neutral_layer,// molecules
		   GridLayer& ion_layer, // ions
		   Reaction React, int ncoll_ion_ntrl[],
		   double icoll_dr[], double icoll_dt[], double ncoll_dr[], double ncoll_dt[])
{
  int i_m;

  double W2, W, IW, v_proj_rt, ivp, cos_theta, sin_theta, cos_beta, sin_beta;
  double cos_phi, sin_phi, psi;

  double vr, vt, vz;
  Vec3d v_rel;

  const double M_m=neutral_layer.mass;
  const double M_i=ion_layer.mass;
  const int dt_ntrl=neutral_layer.dt_subcycling;
  const int dt_ion=ion_layer.dt_subcycling;

  const double dtn_over_dti=(double) neutral_layer.dt_subcycling/ion_layer.dt_subcycling;
  const double dti_over_dtn = 1. / dtn_over_dti;

  double E_th, Emin, Emax, Estep, VolFact;
  int Epoints, Eind;
  E_th = React.Eth;
  Emin = React.Emin;
  Estep = React.Estep;
  Epoints = React.N;
  
  for (int i = 0; i < neutral_layer.r_dim; ++i) 
  {
    // Zylindrical volume scaling according  Verboncoer,JCP 2001, equ: 19 + 20
    if (i != 0) { VolFact = 1. / i; }
    else        { VolFact = 6.;    }

    for (int j = 0; j < neutral_layer.z_dim; ++j) 
    {
      auto& neutral_cell = neutral_layer.get_cell( i , j );
      auto& ion_cell = ion_layer.get_cell( i , j );

      if (!neutral_cell.empty() && !ion_cell.empty()) 
      {
	for (int iid = 0; iid < ion_cell.size(); ++iid) 
	{
	  i_m = (int)(neutral_cell.size() * RAND);
	  if (i_m == neutral_cell.size()) { i_m--; }

	  Particle& ion = ion_cell.particles[iid];
	  Particle& neutral = neutral_cell.particles[i_m];

	  v_rel.r = ion.vr * dtn_over_dti - neutral.vr;  // relative velocity on the neutral scale
	  v_rel.t = ion.vt * dtn_over_dti - neutral.vt;  // !!3.09.01 Check it on scaling 4 ion velocity
	  v_rel.z = ion.vz * dtn_over_dti - neutral.vz;

	  W2 = SQU(v_rel.r) + SQU(v_rel.t) + SQU(v_rel.z);  // before scattering
	  W = sqrt(W2);

	  // Linear fit for Cross-section
	  double S_i = 0.;
	  Eind = (int)((W2 - Emin) / Estep);
	  if (Eind < 0)             { S_i = React.CS[0];       }
	  else if (Eind >= Epoints) { S_i = React.CS[Epoints]; }
	  else {                     
	    S_i = React.CS[Eind] + (React.CS[Eind + 1] - React.CS[Eind]) * ((W2 - Emin) / Estep - Eind);
	  }
	  if (  RAND < neutral_cell.size() * W * S_i * VolFact)  // !! Est' collision!! 1-exp(-n_m*W_0*Si*dt_coll)
	  {
	/*    
	     *engcheck += 0.5*m_e*(SQU(ions[i_el + Next_el].vx) + SQU(ions[i_el + Next_el].vy) + SQU(ions[i_el + Next_el].vz));
	     *engcheck += 0.5*M_m*(SQU(molecules[i_m + Next_m].vx*dti) 
			  + SQU(molecules[i_m+ Next_m].vy*dti) + SQU(molecules[i_m + Next_m].vz*dti)) - 0.5*m_e*M_m/(m_e+M_m)*E_th ;
	     
	     (*momcheck).x += m_e*ions[i_el + Next_el].vx + M_m*molecules[i_m +Next_m].vx*dti;
	     (*momcheck).y += m_e*ions[i_el + Next_el].vy + M_m*molecules[i_m +Next_m].vy*dti;
	     (*momcheck).z += m_e*ions[i_el + Next_el].vz + M_m*molecules[i_m +Next_m].vz*dti;
	   */ 
	    W = MAX(W, 1.e-14);  // A chtob na 0 ne delit'
	    IW = 1. / W;

	    v_proj_rt = sqrt(SQU(v_rel.r) + SQU(v_rel.t));  // proekciya na XY
	    v_proj_rt = MAX(v_proj_rt, 1.e-14);		    // ta zhe fignia
	    ivp = 1. / v_proj_rt;

	    cos_theta = v_rel.z * IW;  // theta - ugol mezhdu v_rel i OZ
	    sin_theta = v_proj_rt * IW;

	    cos_beta = v_rel.r * ivp;  // beta - ugol mezhdu proekciey v_rel na XY i os'yu abciss
	    sin_beta = v_rel.t * ivp;

	    /* Sampling scattering angles in C.M. system */
	    if (RAND < 0.5) { cos_phi = 2. * RAND - 1.; }
	    else            { cos_phi = -1.;            }
	    
	    sin_phi = sqrt(1. - SQU(cos_phi));
	    psi = TWOPI * RAND;

	    vz = W * cos_phi;  //   komponenty posle stolknoveniya v povernutoy coor. sys.
	    vr = W * sin_phi * cos(psi);
	    vt = W * sin_phi * sin(psi);

	    v_rel.r -= vr*cos_beta*cos_theta - vt*sin_beta + vz*cos_beta*sin_theta;  // eto  uzhe minus prirascheniya otn. skor.
	    v_rel.t -= vr*sin_beta*cos_theta + vt*cos_beta + vz*sin_beta*sin_theta;  // a ne otn. skorost'
	    v_rel.z -= vz * cos_theta - vr * sin_theta;

	    ion.vr -= M_m / (M_m + M_i) * v_rel.r * dti_over_dtn;
	    ion.vt -= M_m / (M_m + M_i) * v_rel.t * dti_over_dtn;
	    ion.vz -= M_m / (M_m + M_i) * v_rel.z * dti_over_dtn;
#if !NTRL_CONST
	    neutral.vr += M_i / (M_m + M_i) * v_rel.r;
	    neutral.vt += M_i / (M_m + M_i) * v_rel.t;
	    neutral.vz += M_i / (M_m + M_i) * v_rel.z;
#endif
	    /*
	     *engcheck -= 0.5*m_e*(SQU(ions[i_el + Next_el].vx) + SQU(ions[i_el + Next_el].vy)+ SQU(ions[i_el + Next_el].vz));
	     *engcheck -= 0.5*M_m*(SQU(target_v.x) + SQU(target_v.y)+ SQU(target_v.z));
	     
	     (*momcheck).x -= m_e*ions[i_el + Next_el].vx + M_m*target_v.x;
	     (*momcheck).y -= m_e*ions[i_el + Next_el].vy + M_m*target_v.y;
	     (*momcheck).z -= m_e*ions[i_el + Next_el].vz + M_m*target_v.z;
	    */

// collision diagnostics
#if DIAG_COLL
	  ncoll_ion_ntrl[neutral_layer.cell_number(i,j)]++;
	  diag_collisions( ion, icoll_dr, icoll_dt );
	  diag_collisions( neutral, ncoll_dr, ncoll_dt );
#endif

	  }  //  if (W2_0 >= E_th && RAND < n_m*W_0*Si)
	}  // for ( iid = 0; iid < n_el; iid++)
      }  // if ( n_m && n_el )
    }  // for (j=0; j < nz;  j++)
  }
}

/***********************************************************************************

  Function:     coll_el_ntrl_ioniz(...)
  Action:       neutrals-electrons collisssions with production of ions (e + n -> 2e + i )
		1. We do inelastic collision with neutral in which we lose E_th
		2. We do ellastic collision with electron.

***********************************************************************************/

void coll_el_ntrl_ioniz( GridLayer& neutral_layer,
			 GridLayer& electron_layer,
			 GridLayer& ion_layer, 
			 Reaction React, Vec3d *momcheck,
			double *engcheck, int ncoll_ioniz[], double ecoll_dr[],
			double ecoll_dt[], double ncoll_dr[],
			double ncoll_dt[]) 
{
  const double 	m_e = 1.;
  double 	VolFact;
  double 	W2, W, IW, W2_0, W_0, IW_0, v_proj_rt, ivp;
  double 	cos_theta, sin_theta,cos_beta, sin_beta;
  double 	cos_phi, sin_phi, psi;
  double 	vr, vt, vz;
  Vec3d 	v_rel, delta_v_rel, target_v;
  int 		Eind;

  const double M_n=neutral_layer.mass;
  const int dt_ntrl=neutral_layer.dt_subcycling;
  const int dt_ion=ion_layer.dt_subcycling;

  const double dti = 1. / dt_ion;  // Scaling factor 4 ion velocity
  const double dtn = 1. / dt_ntrl;

  //printf("start coll_el_ntrl_ioniz: \n");
  double E_th = React.Eth;	//printf("E_th=%f\n",E_th);
  double Emin = React.Emin;	//printf("Emin=%f\n",Emin);
  //printf("Emax=%f\n",React.Emax);
  double Estep = React.Estep;	//printf("Estep=%f\n",Estep);
  int Epoints = React.N;	//printf("Epoints%f\n",Epoints);


  for (int i = 0; i < neutral_layer.r_dim; ++i) {
    // Zylindrical volume scaling according  Verboncoer,JCP 2001, equ: 19 + 20
    if (i != 0)
      VolFact = 1. / i;
    else
      VolFact = 6.;

    for (int j = 0; j < neutral_layer.z_dim; ++j) {

      auto& neutral_cell = neutral_layer.get_cell( i , j );
      auto& electron_cell = electron_layer.get_cell( i , j );
      //printf("++cell%i: n_n=%i n_e=%i\n",i_g,n_n,n_e);
      
      if (!neutral_cell.empty() && !electron_cell.empty()) {
	if( neutral_cell.size() > NNcell )   fprintf(stderr,"To much neutrals !!! increase NNcell n_n %zd",neutral_cell.size());   

	int number_of_electrons = electron_cell.size();
	int number_of_neutrals = neutral_cell.size();

for (int i_el = 0; i_el < number_of_electrons && neutral_cell.size(); ++i_el)  // check it out
	{ 
	  int k = (int)(neutral_cell.size() * RAND);
	  if (k == neutral_cell.size()) --k;

	  Particle& electron = electron_cell.particles[i_el];
	  Particle& neutral = neutral_cell.particles[k];

          // relativ velocity (TODO: !!3.09.01 Check it on scaling 4 ion velocity)
	  v_rel.r = electron.vr - neutral.vr * dtn;
	  v_rel.t = electron.vt - neutral.vt * dtn;
	  v_rel.z = electron.vz - neutral.vz * dtn;
	  
	  W2_0 = SQU(v_rel.r) + SQU(v_rel.t) + SQU(v_rel.z);  // before scattering
	  W_0 = sqrt(W2_0);
	  // Linear fit for Cross-section
	  double S_i = 0.;
	  if (W2_0 >= E_th) {
	    Eind = (int)((W2_0 - Emin) / Estep);
	    if (Eind < 0)
	      S_i = React.CS[0];
	    else if (Eind >= Epoints)
	      S_i = React.CS[Epoints];
	    else
	      S_i = React.CS[Eind] + (React.CS[Eind + 1] - React.CS[Eind])*((W2_0 - Emin) / Estep - Eind);
	  }


	  //number_of_neutrals here instead of neutral_cell.size() for compatibility w/ old code
	  //TODO check which is correct
	  if (W2_0 >= E_th &&  RAND < number_of_neutrals * W_0 * S_i * VolFact)  // !! Est' collision!! 1-exp(-n_n*W_0*Si*dt_coll)
	  {
	    
	    *engcheck += 0.5 * m_e * (SQU(electron.vr) +
				      SQU(electron.vt) +
				      SQU(electron.vz));
	    *engcheck += 0.5 * M_n * (SQU(neutral.vr * dtn) +
				      SQU(neutral.vt * dtn) +
				      SQU(neutral.vz * dtn)) - 0.5 * m_e * M_n / (m_e + M_n) * E_th;

	    (*momcheck).r += m_e * electron.vr + M_n * neutral.vr * dtn;
	    (*momcheck).t += m_e * electron.vt + M_n * neutral.vt * dtn;
	    (*momcheck).z += m_e * electron.vz + M_n * neutral.vz * dtn;

	    W_0 = MAX(W_0, 1.e-14);
	    IW_0 = 1. / W_0;

	    W2 = W2_0 - E_th;
	    W = sqrt(W2);  // Absolutnaya velichina otn. skorosti

	    delta_v_rel.r = v_rel.r * (W * IW_0 - 1.);  // After energy loss
	    delta_v_rel.t = v_rel.t * (W * IW_0 - 1.);
	    delta_v_rel.z = v_rel.z * (W * IW_0 - 1.);

	    v_rel.r += delta_v_rel.r;
	    v_rel.t += delta_v_rel.t;
	    v_rel.z += delta_v_rel.z;

	    W = MAX(W, 1.e-14);  // A chtob na 0 ne delit'
	    IW = 1. / W;

	    v_proj_rt = sqrt(SQU(v_rel.r) + SQU(v_rel.t));  // proekciya na XY
	    v_proj_rt = MAX(v_proj_rt, 1.e-14);		    // ta zhe fignia
	    ivp = 1. / v_proj_rt;

	    cos_theta = v_rel.z * IW;  // theta - ugol mezhdu v_rel i OZ
	    sin_theta = v_proj_rt * IW;

	    cos_beta = v_rel.r * ivp;  // beta - ugol mezhdu proekciey v_rel na XY i os'yu abciss
	    sin_beta = v_rel.t * ivp;

	    /* Sampling scattering angles in C.M. system */
	    cos_phi = 2. * RAND - 1.;  //FIXME: just 4 test!! Use = (2.+E-2.*pow((1+E), RAND))/E instead
	    sin_phi = sqrt(1. - SQU(cos_phi));
	    psi = TWOPI * RAND;

	    vz = W * cos_phi;  //   komponenty posle stolknoveniya v povernutoy coor. sys.
	    vr = W * sin_phi * cos(psi);
	    vt = W * sin_phi * sin(psi);

	    v_rel.r -= vr * cos_beta * cos_theta - vt * sin_beta + vz * cos_beta * sin_theta;  // eto  uzhe minus prirascheniya otn. skor.
	    v_rel.t -= vr * sin_beta * cos_theta + vt * cos_beta + vz * sin_beta * sin_theta;  // a ne otn. skorost'
	    v_rel.z -= vz * cos_theta - vr * sin_theta;

	    v_rel.r -= delta_v_rel.r;
	    v_rel.t -= delta_v_rel.t;
	    v_rel.z -= delta_v_rel.z;

	    target_v.r = neutral.vr * dtn;
	    target_v.t = neutral.vt * dtn;
	    target_v.z = neutral.vz * dtn;

	    electron.vr -= M_n / (M_n + m_e) * v_rel.r;
	    electron.vt -= M_n / (M_n + m_e) * v_rel.t;
	    electron.vz -= M_n / (M_n + m_e) * v_rel.z;

	    target_v.r += m_e / (M_n + m_e) * v_rel.r;  // a nado li
	    target_v.t += m_e / (M_n + m_e) * v_rel.t;  //
	    target_v.z += m_e / (M_n + m_e) * v_rel.z;  // Dlia diagnostiki - prisvoy neutralam sho- nibud' nelepoe

	    Particle ion;

	    // Ion is born at position of neutral
	    ion.assign_new_id();
	    ion.r = neutral.r;
	    ion.z = neutral.z;

	    ion.vr = target_v.r * dt_ion;  // @@ scale it back
	    ion.vt = target_v.t * dt_ion;
	    ion.vz = target_v.z * dt_ion;

#if DIAG_COLL
	    ion.coll_r = 0.;
	    ion.coll_t = nstep;
#endif

	    store_particle_properties( &ion, IONIZATION_ORIGIN );

	    *engcheck -= 0.5 * (M_n - m_e) * (SQU(ion.vr * dti) +
					      SQU(ion.vt * dti) +
					      SQU(ion.vz * dti));

	    (*momcheck).r -= (M_n - m_e) * ion.vr * dti;
	    (*momcheck).t -= (M_n - m_e) * ion.vt * dti;
	    (*momcheck).z -= (M_n - m_e) * ion.vz * dti;

	    // put the new ion into the ion layer
	    ion_layer.store( ion );

	    /******  electron-neutral inelastic collision with a loss of E_th energy and birth of the ion  is done !
	             now we'll do electron-electron elastic collission  ******/

	    v_rel.r = electron.vr - target_v.r;  // snova relative velocity
	    v_rel.t = electron.vt - target_v.t;
	    v_rel.z = electron.vz - target_v.z;

	    W2 = SQU(v_rel.r) + SQU(v_rel.t) + SQU(v_rel.z);
	    W = sqrt(W2);	// Absolutnaya velichina otn. skorosti
	    W = MAX(W, 1.e-14);  // A chtob na 0 ne delit'
	    IW = 1. / W;

	    v_proj_rt = sqrt(SQU(v_rel.r) + SQU(v_rel.t));  // proekciya na XY
	    v_proj_rt = MAX(v_proj_rt, 1.e-14);		    // ta zhe fignia
	    ivp = 1. / v_proj_rt;

	    cos_theta = v_rel.z * IW;  // theta - ugol mezhdu v_rel i OZ
	    sin_theta = v_proj_rt * IW;

	    cos_beta = v_rel.r * ivp;  // beta - ugol mezhdu proekciey v_rel na XY i os'yu abciss
	    sin_beta = v_rel.t * ivp;

	    cos_phi = 2. * RAND - 1.;
	    sin_phi = sqrt(1 - SQU(cos_phi));
	    psi = TWOPI * RAND;

	    vz = W * cos_phi;  // komponenty posle stolknoveniya v povernutoy coor. sys.
	    vr = W * sin_phi * cos(psi);
	    vt = W * sin_phi * sin(psi);

	    v_rel.r -= vr * cos_beta * cos_theta - vt * sin_beta + vz * cos_beta * sin_theta;  // eto  uzhe minus prirascheniya otn. skor.
	    v_rel.t -= vr * sin_beta * cos_theta + vt * cos_beta + vz * sin_beta * sin_theta;  // a ne otn. skorost'
	    v_rel.z -= vz * cos_theta - vr * sin_theta;

	    electron.vr -= 0.5 * v_rel.r;  // Note inversed sign due inversed delta be4
	    electron.vt -= 0.5 * v_rel.t;  // Rastiapa!!!
	    electron.vz -= 0.5 * v_rel.z;

	    target_v.r += 0.5 * v_rel.r;
	    target_v.t += 0.5 * v_rel.t;
	    target_v.z += 0.5 * v_rel.z;

	    // if( *ne + inzd +1 > NE )   fprintf(stderr,"Slishkom mnogo
	    // ponastalkivali electronov");
	    
	    Particle new_electron;

	    new_electron.assign_new_id();
	    new_electron.r = neutral.r;
	    new_electron.z = neutral.z;
	    new_electron.vr = target_v.r;
	    new_electron.vt = target_v.t;
	    new_electron.vz = target_v.z;
#if DIAG_COLL
	    new_electron.coll_r = 0.;
	    new_electron.coll_t = nstep;
#endif

	    *engcheck -= 0.5 * m_e *(SQU(electron.vr) +SQU(electron.vt) +SQU(electron.vz));
	    *engcheck -= 0.5 * m_e *(SQU(new_electron.vr) +SQU(new_electron.vt) +SQU(new_electron.vz));

	    (*momcheck).r -= m_e * (electron.vr + new_electron.vr);
	    (*momcheck).t -= m_e * (electron.vt + new_electron.vt);
	    (*momcheck).z -= m_e * (electron.vz + new_electron.vz);

	    // collision diagnostics 
#if DIAG_COLL
	    ncoll_ioniz[i*electron_layer.z_dim+j]++;
	    electron_cell.diagnostic_arrays.coll_ioniz++;
	    diag_collisions( electron, ecoll_dr, ecoll_dt ); 
	    diag_collisions( neutral, ncoll_dr, ncoll_dt ); // for 2nd order ionisation ncoll_ -> icoll_
#endif
	    electron_layer.store(new_electron);

#if !NTRL_CONST
	    unsigned int id = k;
	    neutral_cell.remove( id ) ;
#endif

          }  //  if (W2_0 >= E_th && RAND < n_n*W_0*Si)
	}  // for ( i_el = 0; i_el < n_el; i_el++)
      }  // if ( n_n && n_el )
    }  // for (j=0; j < nz;  j++)
  }  // for (i=0; i < nr;  i++)
}

void calculate_pairs_per_cell(int ordcount_n[], Reaction react, int nr, int nz,
			      long pairs_per_cell[]) {
  double E_th = react.Eth;  // it's vsigma_max now
  double VolFact;
  long total_collisions = 0;

  for (int i = 0; i < nr; ++i) {
    // Zylindrical volume scaling according  Verboncoer,JCP 2001, equ: 19 + 20
    if (i != 0)
      VolFact = 1. / i;
    else
      VolFact = 6.;

    for (int j = 0; j < nz; ++j) {
      int i_g = i * nz + j;
#if 0
      if ( j < 50 ) {
	pairs_per_cell[i*layer.z_dim+j] = 0;	
	continue;
      }
#endif
      int n_n = ordcount_n[i_g];
      long nPairs = 0;

      if (n_n > 1) {
	int n2coll = (int)(E_th * VolFact * n_n * n_n);
	nPairs = n2coll / 2;

	if (nPairs * 2 == n2coll - 1)
	  if (RAND < 0.5) ++nPairs;
      }
      pairs_per_cell[i_g] = nPairs;
      total_collisions += nPairs;
    }
  }

  int mpi_rank = get_rank();
#if USE_MPI
  printf("rank %d total collisions %ld\n", mpi_rank, total_collisions);
#endif

}

/*******************************************************************************/
/* nice features for debugging collisiosns                                     */
void debug_coll_partners(Particle& pa1, Particle& pa2, int mpi_rank){
#if DEBUG_PRTL_COLL
  int r1 =  (int)pa1.r;  int z1 =  (int)pa1.z;
  int r2 =  (int)pa2.r;  int z2 =  (int)pa2.z;

  assert( r1 == r2 && z1==z2 && " colliding particles are not in the same cell" );

  // Check wether cell of particle belongs to current node
  int stop = 0;
  if ( mpi_rank != lookup( z1, r1 ) ) {
    printf("working with particle of other processor\n" );
    printf("rank %d should be rank %d -- r %d z %d -- %d %d \n", mpi_rank, lookup(z1,r1) , i , j, r1, z1 );
    stop = 1;
  }
  if ( mpi_rank != lookup( z2, r2 ) ) {
    printf("working with particle of other processor\n" );
    printf("rank %d should be rank %d -- r %d z %d -- %d %d \n", mpi_rank, lookup(z2,r2) , i , j, r2, z2 );
    stop = 1;
  }
  if ( stop ) {
    exit(-1);
  }

#endif
}

/*******************************************************************************/
/* Arrange the N elements of ARRAY in random order                             */
//void shuffle(Particle *array, int n){
void shuffle(Cell& cell){
  for (int i = cell.size() - 1; i > 0; i--) {
    int j =  RAND*(i+1-std::numeric_limits<double>::min());
    std::swap( cell.particles[j], cell.particles[i] );
  }
}

/*******************************************************************************/
/* ensures particles are collide only once                                     */
void select_unused_particles( int* was_used, int size, int pair_index, int* pa1, int* pa2 ) {
  int id1 = -1;
  int id2 = -1;
  
  if(size <= (2*pair_index)+1 ) { printf(" numb of particles in cell =%i  pair_index=%i\n",size,pair_index); }
  assert( (2*pair_index)+1 < size  && "number of collisions > pairs of particles" );
  id1 = was_used[2*pair_index];
  id2 = was_used[(2*pair_index)+1];

  *pa1 = id1;
  *pa2 = id2;
}

/*******************************************************************************/
/* diagnose collision partners                                                 */
void diagnose_coll_partner( Particle pa[], int pair_index, int ip1, int ip2, 
			    int Next_n, int r, int z,Vec3d v_rel ) 
{
#if DIAG_NN_COLL_PARTNER
  int mpi_rank = get_rank();
  char rank[9];
  snprintf(rank, sizeof(rank), "%02i", mpi_rank);
  FILE *file;
  char Fn[200]; 
  strcpy(Fn, "out/distance_n-n-CollPartners_");
  strcpy(Fn + strlen(Fn), rank);
  strcpy(Fn + strlen(Fn), ".dat");
  file = fopen(Fn, "a");
  assert( file );
  
  //printf(" p_ind=%i:   collide  ip1=%i  with  ip2=%i\n",pair_index,ip1,ip2);

  int index1 = ip1+Next_n;
  int index2 = ip2+Next_n;
  double distance, v_rel;
  if( r<coll_area_r_max && coll_area_r_min<=r && 
      z<coll_area_z_max && coll_area_z_min<=z && pair_index<numb_pairs )
  {
    distance = sqrt( pow(pa[index1].r-pa[index1].r,2) +pow(pa[index1].z-pa[index1].z,2) );
    fprintf(file,"%i %f %f %f %f %f  %i %f %f %f %f %f  %f  %f %f %f\n",
		   ip1, pa[index1].r, pa[index1].z, pa[index1].vr, pa[index1].vz, pa[index1].vt,
		   ip2, pa[index2].r, pa[index2].z, pa[index2].vr, pa[index2].vz, pa[index2].vt,
		   distance, v_rel.r, v_rel.z, v_rel.t );
  }else if( pair_index >= numb_pairs ){
    //printf("diagnostic of n-ncollisons partners is limited to the first %i pairs/cell\n",numb_pairs);
  }
  fclose(file);
#endif
}



/************************************************************************************************
  Neutral - neutral collision
  We apply Null -collision
  !!  E_th = MAX(W*sigma here)
*************************************************************************************************/
void coll_ntrl_ntrl( GridLayer& layer,
		    Reaction React, Vec3d *momcheck, double *engcheck,
		    int ncoll_ntrl_ntrl[], double ncoll_dr[], double ncoll_dt[])

{

  int mpi_rank = get_rank();

  int p_ind, i_n1, i_n2;

  double W2, W, IW, v_proj_rt, ivp, cos_theta, sin_theta, cos_beta, sin_beta;
  double cos_phi, sin_phi, psi;

  double vr, vt, vz;
  Vec3d v_rel;

  double E_th, Emin, Emax, Estep, VolFact;
  int Epoints, Eind;

  E_th = React.Eth;  // it's vsigma_max now
  Emin = React.Emin;
  Emax = React.Emax;
  Estep = React.Estep;
  Epoints = React.N;

  long n2coll_total = 0;
  long coll_total = 0;
  // for statistic diagnostics
  int nCells_LargeNColl = 0; // gives number of cells where n2col > n_n
  float av_n2coll_over_nn = 0.; // gives average of n2coll/n_n

  for (int i = 0; i < layer.r_dim; ++i) {
    // Zylindrical volume scaling according  Verboncoer,JCP 2001, equ: 19 + 20
    if (i != 0)
      VolFact = 1. / i;
    else
      VolFact = 6.;

    for (int j = 0; j < layer.z_dim; ++j) {

      auto& cell = layer.get_cell( i, j );

      if (cell.size() > 1)  // est' s kem stalkivat'sya
      {
        long n2coll = E_th * VolFact * cell.size() * cell.size();
        
	if ( n2coll > cell.size()-1 && j < red_coll_z )
        {
	  //fprintf(stderr," @@ Too big time step by neutral -neutral collisions. Exiting!!!  
	  //r=%d z=%d  n_n=%d  n2coll=%d E_th= %f\n", i, j, n_n,n2coll,E_th );
	  //exit(1);
	  //printf(" @@ Too big time step by neutral -neutral collisions. Exiting!!!  r=%d z=%d  n_n=%lu  n2coll=%lu -> \tn2coll/n_n=%lu\n",i,j,cell.size(),n2coll,n2coll/cell.size());
	  nCells_LargeNColl++;
	  av_n2coll_over_nn += (float)(n2coll)/(float)(cell.size());
#if REDUCE_NTRL_COLL
	  n2coll = cell.size()-1;  
	  //printf("@@ cell (i=%i,j=%i): is set artificially to n_n-1 = %li  (n_n=%i) \n",i,j,n2coll,cell.size());
#endif
	}
	
        
        long nPairs = n2coll / 2;
	if (nPairs * 2 == n2coll - 1){
	  if (RAND < 0.5) { ++nPairs; }
	}
	n2coll_total += nPairs;

#if USE_UNARY_COLLISIONS
	shuffle(cell); 
#endif

	//printf(" cell(%i;%i) -> n_n=%lu   nPairs=%li n2coll=%li\n",i,j,cell.size(),nPairs,n2coll);
	for (p_ind = 0; p_ind < nPairs; ++p_ind) 
	{

#if USE_UNARY_COLLISIONS
	  // each particle collide only once (to ensure random collisons shuffle is done bevore)
	  i_n1 = 2*p_ind;
	  i_n2 = i_n1+1;
#else
	  get_distinct_particle_ids( i_n1, i_n2, cell.size() );
#endif
	  //printf("%i coll partner: in_1=%i in_2=%i (nPairs=%li n_n=%lu)\n",p_ind,i_n1,i_n2,nPairs,cell.size());
	  
	  debug_coll_partners(cell.particles[i_n1], cell.particles[i_n2], mpi_rank);

	  Particle& np1 = cell.particles[i_n1];
	  Particle& np2 = cell.particles[i_n2];

	  // relative velocityi     !!3.09.01 Check it on scaling 4 ion velocity
	  v_rel.r = np1.vr - np2.vr;
	  v_rel.t = np1.vt - np2.vt; 
	  v_rel.z = np1.vz - np2.vz;

	  W2 = SQU(v_rel.r) + SQU(v_rel.t) + SQU(v_rel.z);  // before scattering
	  W = sqrt(W2);

	  // Linear fit for Cross-section
	  double S_i =0.;
	  Eind = (int)((W2 - Emin) / Estep);
	  if (Eind < 0)
	    S_i = React.CS[0];
	  else if (Eind >= Epoints)
	    S_i = React.CS[Epoints];
	  else
	    S_i = React.CS[Eind] +(React.CS[Eind + 1] - React.CS[Eind])*((W2 - Emin) / Estep - Eind);


	  if (RAND < W * S_i / E_th)  // !! Est' collision!! 1-exp(-n_n*W_0*Si*dt_coll)
	  {
	    coll_total++;
	    *engcheck += SQU(np1.vr) +
			 SQU(np1.vt) +
			 SQU(np1.vz);
	    *engcheck += SQU(np2.vr) +
			 SQU(np2.vt) +
			 SQU(np2.vz);

	    (*momcheck).r += np1.vr + np2.vr;
	    (*momcheck).t += np1.vt + np2.vt;
	    (*momcheck).z += np1.vz + np2.vz;
#if 0
	    //if(i==62 && j==0){
	    if(fabs(neutrals[i_n1+Next_n].vr)<0.05 || fabs(neutrals[i_n1+Next_n].vt)<0.05 || fabs(neutrals[i_n1+Next_n].vz)<0.05 || 
		fabs(neutrals[i_n2+Next_n].vr)<0.05 || fabs(neutrals[i_n2+Next_n].vt)<0.05 || fabs(neutrals[i_n2+Next_n].vz)<0.05 )
	    {
	      printf("cell %i %i rank =%i \t i_n1=%i v0=(%.1e ,%.1e ,%.1e) \t i_n2=%i v0=(%.1e ,%.1e ,%.1e)\t v_rel0=(%.1e ,%.1e ,%.1e)\n",
		  i,j,mpi_rank,
		  i_n1,neutrals[i_n1+Next_n].vr, neutrals[i_n1+Next_n].vz, neutrals[i_n1+Next_n].vt, 
		  i_n2,neutrals[i_n2+Next_n].vr, neutrals[i_n2+Next_n].vz, neutrals[i_n2+Next_n].vt, 
		  v_rel.r,v_rel.z,v_rel.r);
	    }
	    //}
#endif
	    W = MAX(W, 1.e-14);  // A chtob na 0 ne delit'
	    IW = 1. / W;

	    v_proj_rt = sqrt(SQU(v_rel.r) + SQU(v_rel.t));  // proekciya na XY
	    v_proj_rt = MAX(v_proj_rt, 1.e-14);		    // ta zhe fignia
	    ivp = 1. / v_proj_rt;

	    cos_theta = v_rel.z * IW;  // theta - ugol mezhdu v_rel i OZ
	    sin_theta = v_proj_rt * IW;

	    cos_beta = v_rel.r * ivp; //beta - ugol mezhdu proekciey v_rel na XY i os'yu abciss
	    sin_beta = v_rel.t * ivp;

	    /* Sampling scattering angles in C.M. system */
	    /* if (RAND < 0.5)
	            cos_phi = 2.*RAND -1. ;
	       else
	            cos_phi = -1.;
	    previous version from konstantin, remove half backward scattering */
	    cos_phi = 2. * RAND - 1;

	    sin_phi = sqrt(1. - SQU(cos_phi));
	    psi = TWOPI * RAND;

	    vz = W * cos_phi;  //   komponenty posle stolknoveniya v povernutoy
			       //   coor. sys.
	    vr = W * sin_phi * cos(psi);
	    vt = W * sin_phi * sin(psi);

	    v_rel.r -= vr * cos_beta * cos_theta - vt * sin_beta + vz * cos_beta *
		    sin_theta;  // eto  uzhe minus prirascheniya otn. skor.
	    v_rel.t -= vr * sin_beta * cos_theta + vt * cos_beta +
		       vz * sin_beta * sin_theta;  // a ne otn. skorost'
	    v_rel.z -= vz * cos_theta - vr * sin_theta;

	    /*        
	    target_v.x = neutrals[i_n2 + Next_n].vx*dti;   // remove it after debuging
	    target_v.y = neutrals[i_n2 + Next_n].vy*dti;   // remove it after debuging
	    target_v.z = neutrals[i_n2 + Next_n].vz*dti;   // remove it after debuging  
	    */
	    
	    // FIXME diagnose_coll_partner(neutrals, p_ind, i_n1, i_n2, Next_n, i, j, v_rel); 

	    np1.vr -= 0.5 * v_rel.r;
	    np1.vt -= 0.5 * v_rel.t;
	    np1.vz -= 0.5 * v_rel.z;

	    np2.vr += 0.5 * v_rel.r;
	    np2.vt += 0.5 * v_rel.t;
	    np2.vz += 0.5 * v_rel.z;
	

	    
	    *engcheck -= SQU(np1.vr) +
			 SQU(np1.vt) +
			 SQU(np1.vz);
	    *engcheck -= SQU(np2.vr) +
			 SQU(np2.vt) +
			 SQU(np2.vz);

	    (*momcheck).r -= np1.vr + np2.vr;
	    (*momcheck).t -= np1.vt + np2.vt;
	    (*momcheck).z -= np1.vz + np2.vz;

// collision diagnostics
#if DIAG_COLL
	    ncoll_ntrl_ntrl[layer.cell_number(i,j)]++;
	    diag_collisions( np1, ncoll_dr, ncoll_dt );
	    diag_collisions( np2, ncoll_dr, ncoll_dt );
#endif


	  }  //  if (W2_0 >= E_th && RAND < cell.size()*W_0*Si)
	}  // for ( i_n1 = 0; i_n1 < n_el; i_n1++)
      }  // if ( cell.size() > 1 )
    }  // for (j=0; j < nz;  j++)
  }

  
  if (nstep%print_statistics_interval==0){ 
    if(nCells_LargeNColl>0){
      printf("%i Cells with #collisions > #particles  \taverage n2coll/n_n = %f \n",nCells_LargeNColl,av_n2coll_over_nn/nCells_LargeNColl);
    }
    printf("rank %d: %ld collisions (tried collisons=%ld) \n", mpi_rank, coll_total, n2coll_total); 
  }
}

//split the heavy particle pa_start -> pa + pa_leftover
void split_particle( Particle pa_start, double weight_aimed, Particle* pa, Particle* pa_leftover)
{
#if USE_PARTICLE_WEIGHTING
  *pa = pa_start;
  pa->w = weight_aimed;
  *pa_leftover = pa_start;
  pa_leftover->w = pa_start.w - weight_aimed;
  //printf("split pa1 (w=%f) -> pa1 (w=%f) + pa1_leftover (w=%f)\n", pa_start.w, pa->w, pa_leftover->w );
  //printf("pa1: v = (%f,%f,%f) w=%f\n",pa->p.vr,pa->p.vz,pa->p.vt, pa->w);
  //printf("pa1_leftover: v = (%f,%f,%f) w=%f\n",pa_leftover->p.vr,pa_leftover->p.vz,pa_leftover->p.vt, pa_leftover->w);
#endif
}
//merged pa1 + pa1_leftover -> pa1, distribute error in energy to pa1
void merge_particle( Particle* pa1, Particle pa1_leftover, Particle pa2, Particle pa2_start, double w1_over_w2)
{
#if USE_PARTICLE_WEIGHTING
  //printf("now merging   w1_over_w2=%f\n",w1_over_w2);
  //printf("pa1: v = (%f,%f,%f) w=%f\n",pa1->p.vr,pa1->p.vz,pa1->p.vt, pa1->w);
  //printf("pa1_leftover: v = (%f,%f,%f) w=%f\n",pa1_leftover.vr,pa1_leftover.vz,pa1_leftover.vt, pa1_leftover.w);
  double vr = (1-w1_over_w2)*pa1_leftover.vr + w1_over_w2 * pa1->p.vr;
  double vz = (1-w1_over_w2)*pa1_leftover.vz + w1_over_w2 * pa1->p.vz;
  double vt = (1-w1_over_w2)*pa1_leftover.vt + w1_over_w2 * pa1->p.vt;
  pa1->p.vr = vr; 
  pa1->p.vz = vz; 
  pa1->p.vt = vt; 
  pa1->w += pa1_leftover.w;
  printf("-> pa1: v = (%f,%f,%f) w=%f\n",pa1->p.vr,pa1->p.vz,pa1->p.vt, pa1->w);

  // energy error
  double v1_start2 = SQU(pa1_leftover.vr)+SQU(pa1_leftover.vz)+SQU(pa1_leftover.vr);
  double v2_start2 = SQU(pa2_start.vr)+SQU(pa2_start.vz)+SQU(pa2_start.vr);
  double v1_2 = SQU(pa1->p.vr)+SQU(pa1->p.vz)+SQU(pa1->p.vr);
  double v2_2 = SQU(pa2.vr)+SQU(pa2.vz)+SQU(pa2.vr);
  //double energy_error = pa1->w*(v1_start2-v1_2) + pa2_start.w*(v2_start2 - v2_2) ;
  double energy_error = v1_start2 + v2_start2 - v1_2 - v2_2 ;
  // distribute energy to pa1 acooring component sizes of v
  double alpha = sqrt(1 + ( fabs(energy_error) / ( pa1->w * v1_2) ));
  printf("energy error =%f -> \talpha = %f\n",energy_error,alpha);
  pa1->p.vr = vr * alpha; 
  pa1->p.vz = vz * alpha; 
  pa1->p.vt = vt * alpha; 
  printf(".. -> pa1: v = (%f,%f,%f) w=%f\n",pa1->p.vr,pa1->p.vz,pa1->p.vt, pa1->w);
#endif
}

/************************************************************************************************
  Neutral - neutral collision   using different weighted particles
  We apply Null -collision
  !!  E_th = MAX(W*sigma here)
*************************************************************************************************/
void coll_ntrl_ntrl_weighted(Particle neutrals[], int ordcount_n[], int nr, int nz,
		    Reaction React, Vec3d *momcheck, double *engcheck,
		    int ncoll_ntrl_ntrl[], double ncoll_dr[], double ncoll_dt[])

{
#if USE_PARTICLE_WEIGHTING
  // for statistic diagnostics
  int   nCells_LargeNColl = 0; // gives number of cells where n2col > n_n
  float av_n2coll_over_nn = 0.; // gives average of n2coll/n_n

  // reaction parameter
  double E_th = React.Eth;  // it's now (v_r *sigma)_max 
  double Emin = React.Emin;
  double Emax = React.Emax;
  double Estep = React.Estep;
  int    Epoints = React.N;

  int  Next_n = 0;
  long n2coll_total = 0;
  long coll_total = 0;
  
  
  for (int i = 0; i < nr; ++i) 
  {
    double VolFact;
    // Zylindrical volume scaling according  Verboncoer,JCP 2001, equ: 19 + 20
    if (i != 0) { VolFact = 1./i; }
    else        { VolFact = 6.;   } 

    for (int j = 0; j < nz; ++j) 
    {
#if 0
       if ( j < 9 ) continue;
#endif
      int i_g = i * nz + j;
      int n_n = ordcount_n[i_g];
      // number of real particles / N_SP in the cell
      double N_n = 0; 
      for (int k = 0; k < n_n; ++k){ N_n += neutrals[k].w; }
      // minimal weight in the current cell
      double w_min_cell = coll_fac_ntrl_ntrl;
      for (int k = 0; k < n_n; ++k){
	if( neutrals[k].w < w_min_cell) { w_min_cell = neutrals[k].w; }
      }

      if (n_n > 1)  // est' s kem stalkivat'sya
      {

        long n2coll = E_th * VolFact * N_n * (N_n / w_min_cell); // TODO: is that realy right?? 
	// have to be (N_n/w_min_cell)*(N_n/cell volume) *dt_coll*Delta t*(v_r*sigma)_max
	// without weighting was: 
	// long n2coll = E_th * VolFact * n_n * n_n;

	if ( n2coll > n_n-1 )
        {
	  //fprintf(stderr," @@ Too big time step by neutral -neutral collisions. Exiting!!!  
	  //r=%d z=%d  n_n=%d  n2coll=%d E_th= %f\n", i, j, n_n,n2coll,E_th );
	  //exit(1);
	  //printf(" @@ Too big time step by neutral -neutral collisions. Exiting!!!  r=%d z=%d  n_n=%d  n2coll=%d -> \tn2coll/n_n=%d\n",i,j,n_n,n2coll,n2coll/n_n);
	  nCells_LargeNColl++;
	  av_n2coll_over_nn += (float)(n2coll)/(float)(n_n);
#if REDUCE_NTRL_COLL
	  n2coll = n_n-1;  //printf("@@ cell (i=%i,j=%i): is set artificially to n_n-1 = %i  (n_n=%i) \n",i,j,n2coll,n_n);
#endif
	}
	
        long nPairs = n2coll / 2;
	if (nPairs * 2 == n2coll - 1){
	  if (RAND < 0.5) { ++nPairs; }
	}
	n2coll_total += nPairs;

	//printf(" cell(%i;%i) -> n_n=%i   nPairs=%i n2coll=%i\n",i,j,n_n,nPairs,n2coll);
	for (int p_ind = 0; p_ind < nPairs; ++p_ind) 
	{

	  int    i_n1, i_n2;
#if USE_UNARY_COLLISIONS
	  // in the HEMP runs, with strong directed source, order function gives 
	  // ordering of the particles also inside the cells = this is better for DSMC 
	  i_n1 = 2*p_ind;
	  i_n2 = (2*p_ind)+1;
#else
	  get_distinct_particle_ids( i_n1, i_n2, n_n );
#endif
	  //printf("pa1: w=%f i=%i r=%f z=%f \n", neutrals[i_n1 + Next_n].w,i_n1,neutrals[i_n1 + Next_n].r,neutrals[i_n1 + Next_n].z);
	  //printf("pa2: w=%f i=%i r=%f z=%f \n", neutrals[i_n2 + Next_n].w,i_n2,neutrals[i_n2 + Next_n].r,neutrals[i_n2 + Next_n].z);
	  // particle 1 is the light one where 2 is the heavieri, otherwise switch
	  if( neutrals[i_n1 + Next_n].w > neutrals[i_n2 + Next_n].w ){
	    int temp_i = i_n1;
	    i_n1 = i_n2;
	    i_n2 = temp_i;
	  }else if ( neutrals[i_n1 + Next_n].w == neutrals[i_n2 + Next_n].w ){
	    //printf("Both colliding partners have the same weight w=%f\n",neutrals[i_n1 + Next_n].w);
	  }
	  double w_min = neutrals[i_n1 + Next_n].w; 
	  double w1_over_w2 =  neutrals[i_n1 + Next_n].w /  neutrals[i_n2 + Next_n].w; 
          
	  //diagnose_coll_partner(neutrals, p_ind, i_n1, i_n2, Next_n, i, j); 
	 

	  // relative velocityi     !!3.09.01 Check it on scaling 4 ion velocity
	  Vec3d  v_rel;
	  v_rel.r = neutrals[i_n1 + Next_n].vr -neutrals[i_n2 + Next_n].vr;
	  v_rel.t = neutrals[i_n1 + Next_n].vt -neutrals[i_n2 + Next_n].vt; 
	  v_rel.z = neutrals[i_n1 + Next_n].vz -neutrals[i_n2 + Next_n].vz;

	  double W2 = SQU(v_rel.r) + SQU(v_rel.t) + SQU(v_rel.z);  // before scattering
	  double W = sqrt(W2);

	  // Linear fit for Cross-section
	  double S_i =0.;
	  int Eind = (int)((W2 - Emin) / Estep);
	  if (Eind < 0)             { S_i = React.CS[0]; }
	  else if (Eind >= Epoints) { S_i = React.CS[Epoints];}	  
	  else{
	    S_i = React.CS[Eind] +(React.CS[Eind + 1] - React.CS[Eind])*((W2 - Emin) / Estep - Eind);
	  }
	  
	  double collisons_propabillity =  W * (S_i/E_th) * (w_min_cell/w_min);
	  if (RAND < collisons_propabillity )  // !! Est' collision!! 1-exp(-n_n*W_0*Si*dt_coll)
	  {
	    Particle pa2 = neutrals[i_n2 + Next_n];
	    Particle pa1, pa1_leftover;
	    //if( neutrals[i_n1 + Next_n].w == neutrals[i_n2 + Next_n].w){
	    //   pa1 = neutrals[i_n1 + Next_n];
	    //}else{ //split the heavy particle 1
	      split_particle(neutrals[i_n1 + Next_n], neutrals[i_n2 + Next_n].w, &pa1, &pa1_leftover); //TODO 
	    //}

	    coll_total++;
	    *engcheck += SQU(pa1.vr) + SQU(pa1.vt) + SQU(pa1.vz);
	    *engcheck += SQU(pa2.vr) + SQU(pa2.vt) + SQU(pa2.vz);
	    (*momcheck).r += pa1.vr + pa2.vr;
	    (*momcheck).t += pa1.vt + pa2.vt;
	    (*momcheck).z += pa1.vz + pa2.vz;

	    W = MAX(W, 1.e-14);  // A chtob na 0 ne delit'
	    double IW = 1. / W;

	    double v_proj_rt = sqrt(SQU(v_rel.r) + SQU(v_rel.t));  // proekciya na XY
	    v_proj_rt = MAX(v_proj_rt, 1.e-14);		    // ta zhe fignia
	    double ivp = 1. / v_proj_rt;

	    double cos_theta = v_rel.z * IW;  // theta - ugol mezhdu v_rel i OZ
	    double sin_theta = v_proj_rt * IW;

	    double cos_beta = v_rel.r * ivp; //beta - ugol mezhdu proekciey v_rel na XY i os'yu abciss
	    double sin_beta = v_rel.t * ivp;

	    /* Sampling scattering angles in C.M. system */
	    /* if (RAND < 0.5){ cos_phi = 2.*RAND -1.; }
	       else           { cos_phi = -1.; }
	    previous version from konstantin, remove half backward scattering */
	    double cos_phi = 2. * RAND - 1;

	    double sin_phi = sqrt(1. - SQU(cos_phi));
	    double psi = TWOPI * RAND;

	    double vz = W * cos_phi;  //   komponenty posle stolknoveniya v povernutoy coor. sys.
	    double vr = W * sin_phi * cos(psi);
	    double vt = W * sin_phi * sin(psi);

	    v_rel.r -= vr * cos_beta * cos_theta - vt * sin_beta + vz * cos_beta * sin_theta;  // eto  uzhe minus prirascheniya otn. skor.
	    v_rel.t -= vr * sin_beta * cos_theta + vt * cos_beta + vz * sin_beta * sin_theta;  // a ne otn. skorost'
	    v_rel.z -= vz * cos_theta - vr * sin_theta;

	    /* Vec3D target_v;
	    target_v.x = pa2.vx*dti;   // remove it after debuging
	    target_v.y = pa2.vy*dti;   // remove it after debuging
	    target_v.z = pa2.vz*dti;   // remove it after debuging  
	    */
 
	    //diagnose_coll_partner(neutrals, p_ind, i_n1, i_n2, Next_n, i, j, v_rel); 

	    pa1.vr -= 0.5 * v_rel.r;
	    pa1.vt -= 0.5 * v_rel.t;
	    pa1.vz -= 0.5 * v_rel.z;

	    pa2.vr += 0.5 * v_rel.r;
	    pa2.vt += 0.5 * v_rel.t;
	    pa2.vz += 0.5 * v_rel.z;

	    //if(neutrals[i_n1 + Next_n].w != neutrals[i_n2 + Next_n].w)
	    //{
	      //merged pa1 with pa1_leftover -> pa1, distribute error in energy to pa1
	      merge_particle( &pa1, pa1_leftover, pa2, neutrals[i_n2 + Next_n], w1_over_w2); //TODO
	    //}

	    *engcheck -= SQU(pa1.vr) + SQU(pa1.vt) + SQU(pa1.vz);
	    *engcheck -= SQU(pa2.vr) + SQU(pa2.vt) + SQU(pa2.vz);
	    (*momcheck).r -= pa1.vr + pa2.vr;
	    (*momcheck).t -= pa1.vt + pa2.vt;
	    (*momcheck).z -= pa1.vz + pa2.vz;

	    neutrals[i_n1 + Next_n] = pa1;
	    neutrals[i_n2 + Next_n] = pa2;

#if DIAG_COLL
	    // collision diagnostics
	    ncoll_ntrl_ntrl[i_g]++;
	    diag_collisions( neutrals, i_n1, Next_n, ncoll_dr, ncoll_dt );
	    diag_collisions( neutrals, i_n2, Next_n, ncoll_dr, ncoll_dt );
#endif


	  }  //  if (W2_0 >= E_th && RAND < n_n*W_0*Si)
	}  // for (int p_ind = 0; p_ind < nPairs; ++p_ind) 

      }  // if ( n_n > 1 )

      Next_n += n_n;
    }  // for (int j=0; j < nz;  j++)
  } // for (int j=0; i < nr;  i++)

  
  int mpi_rank = get_rank();
  if (nstep%print_statistics_interval==0){ 
    if(nCells_LargeNColl>0){
      printf("%i Cells with #collisions > #particles  \taverage n2coll/n_n = %f \n",nCells_LargeNColl,av_n2coll_over_nn/nCells_LargeNColl);
    }
    printf("rank %d: %ld collisions (tried collisons=%ld) \n", mpi_rank, coll_total, n2coll_total); 
  }
#endif
}


/***********************************************************************************

  Function:     coll_el_ntrl_ioniz2(...)
  Action:       neutrals-electrons collisssions with production of 
                    double charged ions (e + n -> 3e + i )
                1. We do inelastic collision with neutral in which we lose E_th
                2. We do 2 ellastic collision of electrons.
  
***********************************************************************************/
void coll_el_ntrl_ioniz2( GridLayer& neutral_layer, GridLayer& electron_layer, GridLayer& ion_layer,
			  Reaction React, Vec3d *momcheck, double *engcheck, int ncoll_double_ioniz[],
			double ecoll_dr[],double ecoll_dt[], double ncoll_dr[],double ncoll_dt[] )
{
  const double 	m_e =electron_layer.mass;
  const double  M_n = neutral_layer.mass;
  double  	VolFact;
  double  	W2, W, IW, W2_0,  W_0, IW_0, v_proj_rt, ivp, cos_theta, sin_theta, cos_beta, sin_beta;
  double  	cos_phi, sin_phi, psi;
  double 	vr, vt, vz;
  Vec3d 	v_rel, delta_v_rel, target_v;

  //printf("start coll_el_ntrl_ioniz2 \n");
  double E_th  = React.Eth;  	//printf("E_th=%f\n",E_th);
  double Emin  = React.Emin; 	//printf("Emin=%f\n",Emin);
  double Emax  = React.Emax; 	//printf("Emax=%f\n",Emax);
  double Estep = React.Estep;	//printf("Estep=%f\n",Estep);
  int    Epoints = React.N;	//printf("Epoints=%d\n",Epoints);

  double dti  = 1./dt_ion;      // Scaling factor 4 ion velocity
  double dtn  = 1./dt_ntrl;

  for (int i=0; i < neutral_layer.r_dim;  ++i)
  {
   // Zylindrical volume scaling according  Verboncoer,JCP 2001, equ: 19 + 20
   if ( i != 0 )   VolFact = 1./i;  
   else            VolFact = 6.;
     
   for (int j=0; j < neutral_layer.z_dim;  ++j)
   {

     auto& neutral_cell = neutral_layer.get_cell( i , j ) ;
     auto& ion_cell = ion_layer.get_cell( i , j ) ;
     auto& electron_cell = electron_layer.get_cell( i , j ) ;

     if ( !neutral_cell.empty() && !electron_cell.empty() )
     { 
	if( neutral_cell.size() > NNcell )   fprintf(stderr,"To much neutrals !!! increase NNcell n_n %zu", neutral_cell.size() );   
	    
	int number_of_electrons = electron_cell.size();
        for ( int i_el = 0; i_el < number_of_electrons && neutral_cell.size() ; ++i_el)    // check it out
        {
           int k = (int)( neutral_cell.size() *RAND);
           if (k == neutral_cell.size()) --k;              

	   Particle& neutral = neutral_cell.particles[k];
	   Particle& electron = electron_cell.particles[i_el];

          v_rel.r= electron.vr - neutral.vr*dtn;      // relative velocity
          v_rel.t= electron.vt - neutral.vt*dtn;                        // !!3.09.01 Check it on scaling 4 ion velocity
          v_rel.z= electron.vz - neutral.vz*dtn;

          W2_0 = SQU(v_rel.r) + SQU(v_rel.t) + SQU(v_rel.z);                // before scattering
          W_0  = sqrt(W2_0);

	  // Linear fit for Cross-section
          double S_i = 0.; 
	  if ( W2_0 >= E_th )
          {
	    int Eind = (int)((W2_0 - Emin)/Estep);
            if ( Eind < 0)
               S_i = React.CS[0];
            else if ( Eind >= Epoints)
               S_i = React.CS[Epoints];
            else
               S_i = React.CS[Eind] + (React.CS[Eind +1] - React.CS[Eind])*((W2_0 - Emin)/Estep - Eind);
           }
	   //printf("** Energy before scattering W2_0=%f (E_th=%f) -> S_i=%f\n",W2_0,E_th,S_i);


	  // TODO this used n_n before which might be wrong because it was never decreased when neutrals were erased
           if (W2_0 >= E_th && RAND < neutral_cell.size() *W_0*S_i*VolFact)      // !! Est' collision!! 1-exp(-n_n*W_0*Si*dt_coll)
           {
 	      *engcheck += 0.5*m_e * ( SQU(electron.vr) + 
				       SQU(electron.vt) + 
				       SQU(electron.vz) );
 	      *engcheck += 0.5*M_n * ( SQU(neutral.vr * dtn) + 
				       SQU(neutral.vt * dtn) + 
				       SQU(neutral.vz * dtn)) - 0.5*m_e*M_n/(m_e+M_n)*E_th;

  	     (*momcheck).r += m_e*electron.vr + M_n*neutral.vr*dtn;
  	     (*momcheck).t += m_e*electron.vt + M_n*neutral.vt*dtn;
  	     (*momcheck).z += m_e*electron.vz + M_n*neutral.vz*dtn;

             W_0  = MAX( W_0, 1.e-14);
             IW_0 = 1./W_0;

             W2 = W2_0 - E_th;
             W  = sqrt(W2);                                 // Absolutnaya velichina otn. skorosti

             delta_v_rel.r = v_rel.r*(W*IW_0 -1.);                       // After energy loss
             delta_v_rel.t = v_rel.t*(W*IW_0 -1.);
             delta_v_rel.z = v_rel.z*(W*IW_0 -1.);

             v_rel.r += delta_v_rel.r;
             v_rel.t += delta_v_rel.t;
             v_rel.z += delta_v_rel.z;

             W  = MAX( W, 1.e-14);                         // A chtob na 0 ne delit'
             IW = 1./W;

             v_proj_rt = sqrt( SQU(v_rel.r) + SQU(v_rel.t) );     // proekciya na XY
             v_proj_rt = MAX( v_proj_rt, 1.e-14);              // ta zhe fignia
             ivp       = 1./v_proj_rt;

             cos_theta = v_rel.z * IW;               // theta - ugol mezhdu v_rel i OZ
             sin_theta = v_proj_rt*IW;

             cos_beta = v_rel.r*ivp;                // beta - ugol mezhdu proekciey v_rel na XY i os'yu abciss
             sin_beta = v_rel.t*ivp;

             /* Sampling scattering angles in C.M. system */
             cos_phi = 2.*RAND -1. ; //FIXME: just 4 test!! Use = (2.+E-2.*pow((1+E), RAND))/E instead
             sin_phi = sqrt(1. - SQU(cos_phi));
             psi  = TWOPI*RAND;

             vz = W*cos_phi;              //   komponenty posle stolknoveniya v povernutoy coor. sys.
             vr = W*sin_phi*cos(psi);
             vt = W*sin_phi*sin(psi);

             v_rel.r -= vr*cos_beta*cos_theta -vt*sin_beta + vz*cos_beta*sin_theta;    // eto  uzhe minus prirascheniya otn. skor.
             v_rel.t -= vr*sin_beta*cos_theta +vt*cos_beta + vz*sin_beta*sin_theta;    // a ne otn. skorost'
             v_rel.z -= vz*cos_theta-vr*sin_theta;

             v_rel.r -= delta_v_rel.r;
             v_rel.t -= delta_v_rel.t;
             v_rel.z -= delta_v_rel.z;

             target_v.r = neutral.vr*dtn;
             target_v.t = neutral.vt*dtn;
             target_v.z = neutral.vz*dtn;

             electron.vr -= M_n/(M_n + m_e)*v_rel.r;
             electron.vt -= M_n/(M_n + m_e)*v_rel.t;
             electron.vz -= M_n/(M_n + m_e)*v_rel.z;

             target_v.r    += m_e/(M_n + m_e)*v_rel.r;   // a nado li
             target_v.t    += m_e/(M_n + m_e)*v_rel.t;   //
             target_v.z    += m_e/(M_n + m_e)*v_rel.z;    //               Dlia diagnostiki - prisvoy neutralam sho- nibud' nelepoe

	     Particle ion;
             //ions[*ni + inzd_ion].y = neutrals[i_n + Next_n].y;          // sozdali ion
	     ion.assign_new_id();
             ion.r = neutral.r;
             ion.z = neutral.z;

             ion.vr = target_v.r*dt_ion;                 // @@ scale it back
             ion.vt = target_v.t*dt_ion;
             ion.vz = target_v.z*dt_ion;
                 
	     store_particle_properties( &ion, IONIZATION_TWOPLUS_ORIGIN );
	    
	     *engcheck -= 0.5*(M_n-m_e)*(SQU(ion.vr*dti) + SQU(ion.vt*dti) + SQU(ion.vz*dti));

	     (*momcheck).r -= (M_n-m_e)*ion.vr*dti;
	     (*momcheck).t -= (M_n-m_e)*ion.vt*dti;
 	     (*momcheck).z -= (M_n-m_e)*ion.vz*dti;

	     ion_cell.push_back( ion );

	     /******  electron-neutral inelastic collision with a loss of E_th energy
        	and birth of the ion  is done !
     		now we'll do first electron-electron elastic collission                ******/

             v_rel.r= electron.vr -  target_v.r;      // snova relative velocity
             v_rel.t= electron.vt -  target_v.t;
             v_rel.z= electron.vz -  target_v.z;

             W2 = SQU(v_rel.r) + SQU(v_rel.t) + SQU(v_rel.z);
             W  = sqrt(W2);                                  // Absolutnaya velichina otn. skorosti
             W  = MAX( W, 1.e-14);                         // A chtob na 0 ne delit'
             IW = 1./W;

             v_proj_rt = sqrt( SQU(v_rel.r) + SQU(v_rel.t) );     // proekciya na XY
             v_proj_rt = MAX( v_proj_rt, 1.e-14);              // ta zhe fignia
             ivp       = 1./v_proj_rt;

             cos_theta = v_rel.z * IW;               // theta - ugol mezhdu v_rel i OZ
             sin_theta = v_proj_rt*IW;

             cos_beta = v_rel.r*ivp;                // beta - ugol mezhdu proekciey v_rel na XY i os'yu abciss
             sin_beta = v_rel.t*ivp;

             cos_phi = 2.*RAND -1. ;
             sin_phi = sqrt(1 - SQU(cos_phi));
             psi  = TWOPI*RAND;

             vz = W*cos_phi;              //   komponenty posle stolknoveniya v povernutoy coor. sys.
             vr = W*sin_phi*cos(psi);
             vt = W*sin_phi*sin(psi);

             v_rel.r -= vr*cos_beta*cos_theta -vt*sin_beta + vz*cos_beta*sin_theta;    // eto  uzhe minus prirascheniya otn. skor.
             v_rel.t -= vr*sin_beta*cos_theta +vt*cos_beta + vz*sin_beta*sin_theta;    // a ne otn. skorost'
             v_rel.z -= vz*cos_theta-vr*sin_theta;

             electron.vr -= 0.5*v_rel.r;         // Note inversed sign due inversed delta be4
             electron.vt -= 0.5*v_rel.t;         // Rastiapa!!!
             electron.vz -= 0.5*v_rel.z;

             target_v.r    += 0.5*v_rel.r;
             target_v.t    += 0.5*v_rel.t;
             target_v.z    += 0.5*v_rel.z;


	     //if( *ne + inzd +1 > NE )   fprintf(stderr,"Slishkom mnogo ponastalkivali electronov");
                
	     Particle new_electron_1;

	     new_electron_1.assign_new_id();
             new_electron_1.r = neutral.r;
             new_electron_1.z = neutral.z;
             new_electron_1.vr = target_v.r;
             new_electron_1.vt = target_v.t;
             new_electron_1.vz = target_v.z;

	     electron_cell.push_back( new_electron_1 );
             
	     /**** now we'll do second electron-electron elastic collission          ******/

             v_rel.r= electron.vr -  target_v.r;      // snova relative velocity
             v_rel.t= electron.vt -  target_v.t;
             v_rel.z= electron.vz -  target_v.z;

             W2 = SQU(v_rel.r) + SQU(v_rel.t) + SQU(v_rel.z);
             W  = sqrt(W2);                                  // Absolutnaya velichina otn. skorosti
             W  = MAX( W, 1.e-14);                         // A chtob na 0 ne delit'
             IW = 1./W;

             v_proj_rt = sqrt( SQU(v_rel.r) + SQU(v_rel.t) );     // proekciya na XY
             v_proj_rt = MAX( v_proj_rt, 1.e-14);              // ta zhe fignia
             ivp       = 1./v_proj_rt;

             cos_theta = v_rel.z * IW;               // theta - ugol mezhdu v_rel i OZ
             sin_theta = v_proj_rt*IW;

             cos_beta = v_rel.r*ivp;                // beta - ugol mezhdu proekciey v_rel na XY i os'yu abciss
             sin_beta = v_rel.t*ivp;

             cos_phi = 2.*RAND -1. ;
             sin_phi = sqrt(1 - SQU(cos_phi));
             psi  = TWOPI*RAND;

             vz = W*cos_phi;              //   komponenty posle stolknoveniya v povernutoy coor. sys.
             vr = W*sin_phi*cos(psi);
             vt = W*sin_phi*sin(psi);

             v_rel.r -= vr*cos_beta*cos_theta -vt*sin_beta + vz*cos_beta*sin_theta;    // eto  uzhe minus prirascheniya otn. skor.
             v_rel.t -= vr*sin_beta*cos_theta +vt*cos_beta + vz*sin_beta*sin_theta;    // a ne otn. skorost'
             v_rel.z -= vz*cos_theta-vr*sin_theta;

             electron.vr -= 0.5*v_rel.r;         // Note inversed sign due inversed delta be4
             electron.vt -= 0.5*v_rel.t;         // Rastiapa!!!
             electron.vz -= 0.5*v_rel.z;

             target_v.r    += 0.5*v_rel.r;
             target_v.t    += 0.5*v_rel.t;
             target_v.z    += 0.5*v_rel.z;


	     //if( *ne + inzd +1 > NE )   fprintf(stderr,"Slishkom mnogo ponastalkivali electronov");
                
	     Particle new_electron_2;
	     new_electron_2.assign_new_id();
             new_electron_2.r = neutral.r;
             new_electron_2.z = neutral.z;
             new_electron_2.vr = target_v.r;
             new_electron_2.vt = target_v.t;
             new_electron_2.vz = target_v.z;
                 
	     electron_cell.push_back( new_electron_2 );


	     *engcheck -= 0.5*m_e*(SQU(electron.vr) + SQU(electron.vt) + SQU(electron.vz));
	     *engcheck -= 0.5*m_e*(SQU(new_electron_1.vr) + SQU(new_electron_1.vt) + SQU(new_electron_1.vz));
 	     *engcheck -= 0.5*m_e*(SQU(new_electron_2.vr) + SQU(new_electron_2.vt) + SQU(new_electron_2.vz));

	    (*momcheck).r -= m_e * ( electron.vr + 
    				     new_electron_1.vr +
				     new_electron_2.vr );
	    (*momcheck).t -= m_e * ( electron.vt + 
    				     new_electron_1.vt +
			    	     new_electron_2.vt );
	    (*momcheck).z -= m_e * ( electron.vz + 
				     new_electron_1.vz +
				     new_electron_2.vz );

#if !NTRL_CONST 
	    unsigned int id = k;
	    neutral_cell.remove( id );
#endif

            // collision diagnostics 
#if DIAG_COLL
            ncoll_double_ioniz[i*ion_layer.z_dim+j]++;
            diag_collisions( electron, ecoll_dr, ecoll_dt ); 
            diag_collisions( neutral, ncoll_dr, ncoll_dt ); 
#endif
               
	   }   //  if (W2_0 >= E_th && RAND < n_n*W_0*Si)
        }  //for (int  i_el = 0; i_el < n_el; i_el++)

     }  // end: if ( n_n && n_el )
    }   //end: for (int j=0; j < nz;  j++)
  }   //end: for (int i=0; i < nr;  i++)
 }

 
/***********************************************************************************

  Function:     coll_el_ntrl_ioniz(...)
  Action:       neutrals-electrons collisssions with production of ions (e + n -> 2e + i )
                1. We do inelastic collision with neutral in which we lose E_th
                2. We do ellastic collision with electron.
   
                ++ Calcullation of ionization rate density Pion 
***********************************************************************************/
 

 void coll_el_ntrl_ioniz_Pion( Particle  neutrals[],  int *nn,  int ordcount_ntrls[], double M_n, int dt_ntrl,
                          Particle  electrons[], int *ne,  int ordcount_el[],
                          Particle  ions[],      int *ni,  int dt_ion,
						  int nr,   int nz,  Reaction React,   /* ++ */  int Pion[],
                          Vec3d *momcheck, double *engcheck )
{
 const double m_e =1.;
 double  VolFact;
 double  W2, W, IW, W2_0,  W_0, IW_0, v_proj_rt, ivp, cos_theta, sin_theta, cos_beta, sin_beta;
 double  cos_phi, sin_phi, psi;
 double vr, vt, vz;
 Vec3d 	v_rel, delta_v_rel, target_v;
 char   is_there[NNcell];

 double E_th  = React.Eth;
 double Emin  = React.Emin;
 double Emax  = React.Emax;
 int Estep = React.Estep;
 int Epoints = React.N;
  
 int inzd = 0;  
 int Next_el = 0;
 int Next_n  = 0;
 int lost_inzd_ntrl = 0;
  
 double dti  = 1./dt_ion;      // Scaling factor 4 ion velocity
 double dtn  = 1./dt_ntrl;

 for (int i=0; i < nr;  ++i)
 {
   // Zylindrical volume scaling according  Verboncoer,JCP 2001, equ: 19 + 20
   if ( i != 0 )   VolFact = 1./i;  
   else            VolFact = 6.;
     
   for (int j=0; j < nz;  ++j)
   {
     int i_g = i*nz + j;	
	      
     int n_n  = ordcount_ntrls[i_g];
     int n_el = ordcount_el[i_g];

     if ( n_n && n_el )
      { 
	    if( n_n > NNcell )   fprintf(stderr,"To much neutrals !!! increase NNcell n_n %d",n_n);   
	    
        for ( int i_n  = 0; i_n < n_n;  ++i_n)
         {
           is_there[i_n]= 1;
         }
         
        int N2coll = n_n;
        
        for ( int i_el = 0; i_el < n_el && N2coll ; ++i_el)    // check it out
         {
           int k = (int)(N2coll*RAND);
           if (k == N2coll) --k;              

           int i_n = 0;             
           int i_next   = is_there[i_n];

          while ( i_next <= k )             
           {
            i_next += is_there[++i_n];    
           }

          v_rel.r= electrons[i_el + Next_el].vr - neutrals[i_n + Next_n].vr*dtn;      // relative velocity
          v_rel.t= electrons[i_el + Next_el].vt - neutrals[i_n + Next_n].vt*dtn;                        // !!3.09.01 Check it on scaling 4 ion velocity
          v_rel.z= electrons[i_el + Next_el].vz - neutrals[i_n + Next_n].vz*dtn;

          W2_0 = SQU(v_rel.r) + SQU(v_rel.t) + SQU(v_rel.z);                // before scattering
          W_0  = sqrt(W2_0);
	  
	  // Linear fit for Cross-section
	  double S_i = 0.;
	  if ( W2_0 >= E_th )
          {
               int Eind = (int)((W2_0 - Emin)/Estep);
               if ( Eind < 0)
                   S_i = React.CS[0];
               else if ( Eind >= Epoints)
                   S_i = React.CS[Epoints];
               else
                   S_i = React.CS[Eind] + (React.CS[Eind +1] - React.CS[Eind])*((W2_0 - Emin)/Estep - Eind);
              }


             if (W2_0 >= E_th && RAND < n_n*W_0*S_i*VolFact)      // !! Est' collision!! 1-exp(-n_n*W_0*Si*dt_coll)
               {
 *engcheck += 0.5*m_e*(SQU(electrons[i_el + Next_el].vr) + SQU(electrons[i_el + Next_el].vt) + SQU(electrons[i_el + Next_el].vz));
 *engcheck += 0.5*M_n*(SQU(neutrals[i_n + Next_n].vr*dtn) + SQU(neutrals[i_n + Next_n].vt*dtn) + SQU(neutrals[i_n + Next_n].vz*dtn)) - 0.5*m_e*M_n/(m_e+M_n)*E_th;

  (*momcheck).r += m_e*electrons[i_el + Next_el].vr + M_n*neutrals[i_n + Next_n].vr*dtn;
  (*momcheck).t += m_e*electrons[i_el + Next_el].vt + M_n*neutrals[i_n + Next_n].vt*dtn;
  (*momcheck).z += m_e*electrons[i_el + Next_el].vz + M_n*neutrals[i_n + Next_n].vz*dtn;

                W_0  = MAX( W_0, 1.e-14);
                IW_0 = 1./W_0;

                W2 = W2_0 - E_th;
                W  = sqrt(W2);                                 // Absolutnaya velichina otn. skorosti

                delta_v_rel.r = v_rel.r*(W*IW_0 -1.);                       // After energy loss
                delta_v_rel.t = v_rel.t*(W*IW_0 -1.);
                delta_v_rel.z = v_rel.z*(W*IW_0 -1.);

                v_rel.r += delta_v_rel.r;
                v_rel.t += delta_v_rel.t;
                v_rel.z += delta_v_rel.z;

                W  = MAX( W, 1.e-14);                         // A chtob na 0 ne delit'
                IW = 1./W;

                v_proj_rt = sqrt( SQU(v_rel.r) + SQU(v_rel.t) );     // proekciya na XY
                v_proj_rt = MAX( v_proj_rt, 1.e-14);              // ta zhe fignia
                ivp       = 1./v_proj_rt;

                cos_theta = v_rel.z * IW;               // theta - ugol mezhdu v_rel i OZ
                sin_theta = v_proj_rt*IW;

                cos_beta = v_rel.r*ivp;                // beta - ugol mezhdu proekciey v_rel na XY i os'yu abciss
                sin_beta = v_rel.t*ivp;

                 /* Sampling scattering angles in C.M. system */
                 cos_phi = 2.*RAND -1. ; //FIXME: just 4 test!! Use = (2.+E-2.*pow((1+E), RAND))/E instead
                 sin_phi = sqrt(1. - SQU(cos_phi));
                 psi  = TWOPI*RAND;

                 vz = W*cos_phi;              //   komponenty posle stolknoveniya v povernutoy coor. sys.
                 vr = W*sin_phi*cos(psi);
                 vt = W*sin_phi*sin(psi);

                 v_rel.r -= vr*cos_beta*cos_theta -vt*sin_beta + vz*cos_beta*sin_theta;    // eto  uzhe minus prirascheniya otn. skor.
                 v_rel.t -= vr*sin_beta*cos_theta +vt*cos_beta + vz*sin_beta*sin_theta;    // a ne otn. skorost'
                 v_rel.z -= vz*cos_theta-vr*sin_theta;

                 v_rel.r -= delta_v_rel.r;
                 v_rel.t -= delta_v_rel.t;
                 v_rel.z -= delta_v_rel.z;

                   target_v.r = neutrals[i_n + Next_n].vr*dtn;
                   target_v.t = neutrals[i_n + Next_n].vt*dtn;
                   target_v.z = neutrals[i_n + Next_n].vz*dtn;

                 electrons[i_el + Next_el].vr -= M_n/(M_n + m_e)*v_rel.r;
                 electrons[i_el + Next_el].vt -= M_n/(M_n + m_e)*v_rel.t;
                 electrons[i_el + Next_el].vz -= M_n/(M_n + m_e)*v_rel.z;

                 target_v.r    += m_e/(M_n + m_e)*v_rel.r;   // a nado li
                 target_v.t    += m_e/(M_n + m_e)*v_rel.t;   //
                 target_v.z    += m_e/(M_n + m_e)*v_rel.z;    //               Dlia diagnostiki - prisvoy neutralam sho- nibud' nelepoe

 // if( *ni + inzd > NI )   fprintf(stderr,"Slishkom mnogo naionizovali ionov %d\n", *ni + inzd);
                 //ions[*ni + inzd].y = neutrals[i_n + Next_n].y;          // sozdali ion
		 ions[*ni + inzd].assign_new_id();
                 ions[*ni + inzd].r = neutrals[i_n + Next_n].r;
                 ions[*ni + inzd].z = neutrals[i_n + Next_n].z;

                 ions[*ni + inzd].vr = target_v.r*dt_ion;                 // @@ scale it back
                 ions[*ni + inzd].vt = target_v.t*dt_ion;
                 ions[*ni + inzd].vz = target_v.z*dt_ion;
                 

*engcheck -= 0.5*(M_n-m_e)*(SQU(ions[*ni + inzd].vr*dti) + SQU(ions[*ni + inzd].vt*dti) + SQU(ions[*ni + inzd].vz*dti));

(*momcheck).r -= (M_n-m_e)*ions[*ni + inzd].vr*dti;
(*momcheck).t -= (M_n-m_e)*ions[*ni + inzd].vt*dti;
(*momcheck).z -= (M_n-m_e)*ions[*ni + inzd].vz*dti;

 /******  electron-neutral inelastic collision with a loss of E_th energy
           and birth of the ion  is done !
            now we'll do electron-electron elastic collission                ******/

                 v_rel.r= electrons[i_el + Next_el].vr -  target_v.r;      // snova relative velocity
                 v_rel.t= electrons[i_el + Next_el].vt -  target_v.t;
                 v_rel.z= electrons[i_el + Next_el].vz -  target_v.z;

                 W2 = SQU(v_rel.r) + SQU(v_rel.t) + SQU(v_rel.z);
                 W  = sqrt(W2);                                  // Absolutnaya velichina otn. skorosti
                 W  = MAX( W, 1.e-14);                         // A chtob na 0 ne delit'
                 IW = 1./W;

                 v_proj_rt = sqrt( SQU(v_rel.r) + SQU(v_rel.t) );     // proekciya na XY
                 v_proj_rt = MAX( v_proj_rt, 1.e-14);              // ta zhe fignia
                 ivp       = 1./v_proj_rt;

                 cos_theta = v_rel.z * IW;               // theta - ugol mezhdu v_rel i OZ
                 sin_theta = v_proj_rt*IW;

                 cos_beta = v_rel.r*ivp;                // beta - ugol mezhdu proekciey v_rel na XY i os'yu abciss
                 sin_beta = v_rel.t*ivp;

                 cos_phi = 2.*RAND -1. ;
                 sin_phi = sqrt(1 - SQU(cos_phi));
                 psi  = TWOPI*RAND;

                 vz = W*cos_phi;              //   komponenty posle stolknoveniya v povernutoy coor. sys.
                 vr = W*sin_phi*cos(psi);
                 vt = W*sin_phi*sin(psi);

                 v_rel.r -= vr*cos_beta*cos_theta -vt*sin_beta + vz*cos_beta*sin_theta;    // eto  uzhe minus prirascheniya otn. skor.
                 v_rel.t -= vr*sin_beta*cos_theta +vt*cos_beta + vz*sin_beta*sin_theta;    // a ne otn. skorost'
                 v_rel.z -= vz*cos_theta-vr*sin_theta;

                 electrons[i_el + Next_el].vr -= 0.5*v_rel.r;         // Note inversed sign due inversed delta be4
                 electrons[i_el + Next_el].vt -= 0.5*v_rel.t;         // Rastiapa!!!
                 electrons[i_el + Next_el].vz -= 0.5*v_rel.z;

                 target_v.r    += 0.5*v_rel.r;
                 target_v.t    += 0.5*v_rel.t;
                 target_v.z    += 0.5*v_rel.z;


 //if( *ne + inzd +1 > NE )   fprintf(stderr,"Slishkom mnogo ponastalkivali electronov");
                
		 electrons[*ne + inzd].assign_new_id();
                 electrons[*ne + inzd].r = neutrals[i_n + Next_n].r;
                 electrons[*ne + inzd].z = neutrals[i_n + Next_n].z;
                 electrons[*ne + inzd].vr = target_v.r;
                 electrons[*ne + inzd].vt = target_v.t;
                 electrons[*ne + inzd].vz = target_v.z;
                 


 *engcheck -= 0.5*m_e*(SQU(electrons[i_el + Next_el].vr) + SQU(electrons[i_el + Next_el].vt) + SQU(electrons[i_el + Next_el].vz));
 *engcheck -= 0.5*m_e*(SQU(electrons[*ne + inzd].vr) + SQU(electrons[*ne + inzd].vt) + SQU(electrons[*ne + inzd].vz));

(*momcheck).r -= m_e*(electrons[i_el + Next_el].vr + electrons[*ne + inzd].vr);
(*momcheck).t -= m_e*(electrons[i_el + Next_el].vt + electrons[*ne + inzd].vt);
(*momcheck).z -= m_e*(electrons[i_el + Next_el].vz + electrons[*ne + inzd].vz);

                  ++inzd;
                  --N2coll;
                  is_there[i_n]= 0;        

                   // ionization rate 
                  Pion[i_g]++; 
             }             //  if (W2_0 >= E_th && RAND < n_n*W_0*Si)
          }  //for ( i_el = 0; i_el < n_el; i_el++)
#if !NTRL_CONST             
          for (int i_n = 0; i_n < n_n; ++i_n)                
           {
             if ( !is_there[i_n] )  lost_inzd_ntrl++ ;
             else neutrals[Next_n + i_n - lost_inzd_ntrl] = neutrals[Next_n + i_n]; 
           }
#endif

     }                    //if ( n_n && n_el )
#if !NTRL_CONST
     else for (int i_n = 0; i_n < n_n; i_n++) { 
	neutrals[Next_n + i_n - lost_inzd_ntrl] = neutrals[Next_n + i_n];
     }
#endif
   Next_el += n_el;
   Next_n += n_n;
     }                    //for (int j=0; j < nz;  j++)
    }                     //for (int i=0; i < nr;  i++)
#if !NTRL_CONST
   *nn -= lost_inzd_ntrl;
#endif
   *ne += inzd;  
   *ni += inzd;   
  }
  
 void neutral_collisions( 
	      double* engcheck, 
	      Vec3d* momcheck, 
	      double* colN_time, 
	      double* orderN_time, 
	      Reaction reaction )
{
  int mpi_rank = get_rank(); 

  auto begin_time = Benchmark::start();

  iiprintf("neutral_collisions: decompose domain ... (rank=%d)\n",mpi_rank);
  decompose_domain( global_grid.layers[NEUTRALS] );
  iiprintf("neutral_collisions: ... ready. (rank=%d)\n",mpi_rank);

#if USE_PARTICLE_WEIGHTING
  for (int i = 0; i <nn; ++i){
    if( neutrals[i].w <= 0 ){ printf("** pa: w=%f i=%i r=%f z=%f \n", neutrals[i].w, i, neutrals[i].r, neutrals[i].z); }
  }
#endif

#if USE_PARTICLE_WEIGHTING
  for (int i = 0; i <nn; ++i){
    if( neutrals[i].w <= 0 ){ printf("*** pa: w=%f i=%i r=%f z=%f \n", neutrals[i].w, i, neutrals[i].r, neutrals[i].z); }
  }
#endif

  *orderN_time += Benchmark::stop(begin_time); 

  // all processes distinct from 0 have to zero the values for reduction
  if ( mpi_rank != 0 ) {
    *engcheck = 0;
    momcheck->r = 0;
    momcheck->t = 0;
    momcheck->z = 0;
    for (int i = 0; i < NG; ++i){
        ncoll_ntrl_ntrl[i] = 0;
    }
    for (int i = 0; i < 20*NBIN; ++i){
      ncoll_dr[i] = ncoll_dt[i] = 0;
    }
  }
  
  // Xe - Xe elastic collision
#if USE_PARTICLE_WEIGHTING
  coll_ntrl_ntrl_weighted(neutrals, ordcount_n, nr, nz, reaction, momcheck,
			  engcheck,ncoll_ntrl_ntrl, ncoll_dr, ncoll_dt);
#else 
  coll_ntrl_ntrl(global_grid.layers[NEUTRALS], reaction, momcheck, engcheck,ncoll_ntrl_ntrl, ncoll_dr, ncoll_dt);
#endif

  *colN_time += Benchmark::stop(begin_time);

  reduce( ncoll_ntrl_ntrl, NG);

  reduce( ncoll_dr, 20*NBIN);
  reduce( ncoll_dt, 20*NBIN);

  reduce(engcheck);
  reduce(&momcheck->r);
  reduce(&momcheck->t);
  reduce(&momcheck->z);
}




 /***********************************************************************************

  Function:     coll_el_diss_attach(...)
  Action:       neutrals-electrons collisssions with production of ions (e + n2 -> i- + n )
                1. We do inelastic collision with neutral in which we lose E_th
                2. We do ellastic collision with electron.
  
***********************************************************************************/ 

void coll_el_diss_attach( GridLayer& neutral_layer, 
                          GridLayer& electron_layer,
                          GridLayer& ion_layer, 
			 // GridLayer& atom_layer,
			  double M_a,   // TODO replace with atom neutral layer
			  Reaction React,
                          Vec3d *momcheck, double *engcheck, int ncoll_ioniz[],
			  double ecoll_dr[], double ecoll_dt[], double ncoll_dr[], double ncoll_dt[] )
{

 const double m_e =1.;

 double  dti, dtn, VolFact;
 double  W2, W, IW, W2_0,  W_0, IW_0, v_proj_rt, ivp, cos_theta, sin_theta, cos_beta, sin_beta;
 double  cos_phi, sin_phi, psi;

 double vr, vt, vz;

 const double M_n=neutral_layer.mass;
 const double M_i=ion_layer.mass;
 const int dt_ntrl=neutral_layer.dt_subcycling;
 const int dt_ion=ion_layer.dt_subcycling;
 

 Vec3d 	v_rel, delta_v_rel, target_v, v_CM;
 
 double E_th, Emin, Emax, Estep, S_i;
 int    Epoints, Eind;


  E_th  = React.Eth;
  Emin  = React.Emin;
  Emax  = React.Emax;
  Estep = React.Estep;
  Epoints = React.N;
  S_i = 0.;
  
  dti  = 1./dt_ion;      // Scaling factor 4 ion velocity
  
  dtn  = 1./dt_ntrl;

 for (int i=0; i < neutral_layer.r_dim;  ++i)
  {
   if ( i != 0 )   VolFact = 1./i;  
   else            VolFact = 6.;
     
   for (int j=0; j < neutral_layer.z_dim;  ++j)
    {
     //i_g = i*nz + j;	
	auto& neutral_cell = neutral_layer.get_cell( i , j );
        auto& electron_cell = electron_layer.get_cell( i , j );

     if (neutral_cell.size() && electron_cell.size()) {
        if( neutral_cell.size() > NNcell )   fprintf(stderr,"To much neutrals !!! increase NNcell n_n %zd",neutral_cell.size());   
	    
        for (unsigned int i_el = 0; i_el < electron_cell.size() && neutral_cell.size() ; ++i_el)    // check it out
         {
           int k = (int)(neutral_cell.size()*RAND);
           if (k == neutral_cell.size()) --k;              

		
          Particle& electron = electron_cell.particles[i_el];
          Particle& neutral = neutral_cell.particles[k];

          v_rel.r= electron.vr - neutral.vr*dtn;      // relative velocity
          v_rel.t= electron.vt - neutral.vt*dtn;                        // !!3.09.01 Check it on scaling 4 ion velocity
          v_rel.z= electron.vz - neutral.vz*dtn;

          W2_0 = SQU(v_rel.r) + SQU(v_rel.t) + SQU(v_rel.z);                // before scattering
          W_0  = sqrt(W2_0);
 // Linear fit for Cross-section
 	  double S_i=0;
             if ( W2_0 >= E_th )
              {
               Eind = (int)((W2_0 - Emin)/Estep);
               if ( Eind < 0)
                   S_i = React.CS[0];
               else if ( Eind >= Epoints)
                   S_i = React.CS[Epoints];
               else
                   S_i = React.CS[Eind] + (React.CS[Eind +1] - React.CS[Eind])*((W2_0 - Emin)/Estep - Eind);
	      }
              S_i *=  VolFact;

             if (W2_0 >= E_th && RAND < neutral_cell.size()*W_0*S_i)      // !! Est' collision!! 1-exp(-neutral_cell.size()*W_0*Si*dt_coll)
               {
 
		   
 *engcheck += 0.5*m_e*(SQU(electron.vr) + SQU(electron.vt) + SQU(electron.vz));
 *engcheck += 0.5*M_n*(SQU(neutral.vr*dtn) + SQU(neutral.vt*dtn) + SQU(neutral.vz*dtn)) - 0.5*m_e*M_n/(m_e+M_n)*E_th;

  (*momcheck).r += m_e*electron.vr + M_n*neutral.vr*dtn;
  (*momcheck).t += m_e*electron.vt + M_n*neutral.vt*dtn;
  (*momcheck).z += m_e*electron.vz + M_n*neutral.vz*dtn;

      v_CM.r += (m_e*electron.vr + M_n*neutral.vr*dtn)/(M_n + m_e);
      v_CM.t += (m_e*electron.vt + M_n*neutral.vt*dtn)/(M_n+ m_e);
      v_CM.z += (m_e*electron.vz + M_n*neutral.vz*dtn)/(M_n + m_e);


                W2 = (W2_0 - E_th)*m_e*M_n*(M_i+M_a)/((M_n + m_e)*M_i*M_a);
                W  = sqrt(W2);                                 // Absolutnaya velichina otn. skorosti O i O-

	     	cos_phi = 2.*RAND -1. ; // just 4 test!! Use = (2.+E-2.*pow((1+E), RAND))/E instead
                sin_phi = sqrt(1. - SQU(cos_phi));
                psi  = TWOPI*RAND;


		v_rel.r = W*sin_phi*cos(psi);
                v_rel.t = W*sin_phi*sin(psi);
                v_rel.z = W*cos_phi;

		Particle ion;
		//Negative Ion is born at position of neutral
		//TODO pm, guba: check if right
		ion.assign_new_id();
		ion.r= neutral.r;
		ion.z= neutral.z;

		ion.vr = (v_CM.r + M_a/(M_a + M_i)*v_rel.r)*dti;
		ion.vt = (v_CM.t + M_a/(M_a + M_i)*v_rel.t)*dti;
		ion.vz = (v_CM.z + M_a/(M_a + M_i)*v_rel.z)*dti;

#if DIAG_COLL
		ion.coll_r = 0.;
		ion.coll_t = nstep;    
#endif
		store_particle_properties(&ion, IONIZATION_ORIGIN);
		ion_layer.store( ion ); 

	//	Particle& atom;
	//	printf("negative ion produced at nr %f nz %f\n",ions[*ni + inzd].r,ions[*ni + inzd].r);
#if !NTRL_CONST		
	/*	atom.r= electron.r;
		atom.z= electron.z;

		atom.vr = (v_CM.r - M_i/(M_a + M_i)*v_rel.r)*dti;
		atom.vt = (v_CM.t - M_i/(M_a + M_i)*v_rel.t)*dti;
		atom.vz = (v_CM.z - M_i/(M_a + M_i)*v_rel.z)*dti;

#if DIAG_COLL
		atom.coll_r = 0.;
		atom.coll_t = nstep;*/
#endif
		//store_particle_properties(&atom, IONIZATION_ORIGIN); TODO check if this is needed here


	*engcheck -= 0.5*(M_i)*(SQU(ion.vr*dti) + SQU(ion.vt*dti) + SQU(ion.vz*dti));
//	*engcheck -= 0.5*(M_a)*(SQU(atom.vr*dti) + SQU(atom.vt*dti) + SQU(atom.vz*dti));

//	(*momcheck).r -= (M_i*ion.vr + M_a*atom.vr)*dti;
//	(*momcheck).t -= (M_i*ion.vt + M_a*atom.vt)*dti;
//	(*momcheck).z -= (M_i*ion.vz + M_a*atom.vz)*dti;

	//put the new ion into the ion layer
 
	 //collision diagnostics
#if DEBUG_COLL
#if DIAG_COLL
	ncoll_ioniz[i*electron_layer.z_dim+j]++;
	diag_collisions(electron, ecoll_dr, ecoll_dt);
	diag_collisions(neutral, ncoll_dr, ncoll_dt);
#endif
#endif

#if !NTRL_CONST
	unsigned int id=k;
	neutral_cell.remove(id);
#endif
	electron_cell.remove(i_el);
	    
	
            }   //  if (W2_0 >= E_th && RAND < neutral_cell.size()*W_0*Si)
          }  //for ( i_el = 0; i_el < n_el; i_el++)        
        }                    //if ( neutral_cell.size() && n_el )
      }                    //for (j=0; j < nz;  j++)
    }                     //for (i=0; i < nr;  i++)
  }
    
/***********************************************************************************

  Function:     coll_el_detachment(...)
  
  Action:       ion-electrons collisssions with production of ions (e + O- -> 2e + O )
                1. We do inelastic collision with neutral in which we lose E_th
                2. We do ellastic collision with electron.
COmment: neutrals=negative ions and ions are O 
  
***********************************************************************************/
 

 void coll_el_detachment( GridLayer& neutral_layer,
                          GridLayer&  electron_layer,
                          GridLayer&  ion_layer,
			  Reaction React,
                          Vec3d *momcheck, double *engcheck, int ncoll_ioniz[],
			  double ecoll_dr[], double ecoll_dt[], double ncoll_dr[], double ncoll_dt[] )
{

 const double m_e =1.;

 double  dti, dtn, VolFact;
 double  W2, W, IW, W2_0,  W_0, IW_0, v_proj_rt, ivp, cos_theta, sin_theta, cos_beta, sin_beta;
 double  cos_phi, sin_phi, psi;

 double vr, vt, vz;
 
 const double M_n=neutral_layer.mass;
 const int dt_ntrl=neutral_layer.dt_subcycling;

 Vec3d 	v_rel, delta_v_rel, target_v;
 
 double E_th, Emin, Emax, Estep, S_i;
 int    Epoints, Eind;


  E_th  = React.Eth;
  Emin  = React.Emin;
  Emax  = React.Emax;
  Estep = React.Estep;
  Epoints = React.N;
  S_i = 0.;

  dti  = 1./dt_ion;      // Scaling factor 4 ion velocity
  dtn  = 1./dt_ntrl;	 // Scaling favtor 4 velocity of neutrals

 for (int i=0; i < neutral_layer.r_dim;  ++i)
  {
   if ( i != 0 )   VolFact = 1./i;  
   else            VolFact = 6.;
     
   for (int j=0; j < neutral_layer.z_dim;  ++j)
    {
     //i_g = i*nz + j;	
      
     auto& neutral_cell=neutral_layer.get_cell(i,j);
     auto& electron_cell=electron_layer.get_cell(i,j);

     if ( neutral_cell.size() && electron_cell.size() )
      { 
	    if( neutral_cell.size() > NNcell )   fprintf(stderr,"To much neutrals !!! increase NNcell");   
	    
       int number_of_electrons = electron_cell.size();
        for (unsigned int i_el = 0; i_el < number_of_electrons && neutral_cell.size() ; ++i_el)    // check it out
         {
           int k = (int)(neutral_cell.size()*RAND);
           if (k == neutral_cell.size()) --k;              

	Particle& electron = electron_cell.particles[i_el];
	Particle& neutral = neutral_cell.particles[k];
     // printf("velocities e r %f e t %f e z %f\n", electrons[i_el + Next_el].vr, electrons[i_el + Next_el].vt, electrons[i_el + Next_el].vz );
          v_rel.r= electron.vr - neutral.vr*dtn;      // relative velocity
          v_rel.t= electron.vt - neutral.vt*dtn;                        // !!3.09.01 Check it on scaling 4 ion velocity
          v_rel.z= electron.vz - neutral.vz*dtn;

          W2_0 = SQU(v_rel.r) + SQU(v_rel.t) + SQU(v_rel.z);                // before scattering
          W_0  = sqrt(W2_0);
 // Linear fit for Cross-section
             if ( W2_0 >= E_th )
              {
               Eind = (int)((W2_0 - Emin)/Estep);
               if ( Eind < 0)
                   S_i = React.CS[0];
               else if ( Eind >= Epoints)
                   S_i = React.CS[Epoints];
               else
                   S_i = React.CS[Eind] + (React.CS[Eind +1] - React.CS[Eind])*((W2_0 - Emin)/Estep - Eind);
              }

              S_i *=  VolFact;

             if (W2_0 >= E_th && RAND < neutral_cell.size()*W_0*S_i)      // !! Est' collision!! 1-exp(-neutral_cell.size()*W_0*Si*dt_coll)
               {
 
		   
 *engcheck += 0.5*m_e*(SQU(electron.vr) + SQU(electron.vt) + SQU(electron.vz));
 *engcheck += 0.5*M_n*(SQU(neutral.vr*dtn) + SQU(neutral.vt*dtn) + SQU(neutral.vz*dtn)) - 0.5*m_e*M_n/(m_e+M_n)*E_th;

  (*momcheck).r += m_e*electron.vr + M_n*neutral.vr*dtn;
  (*momcheck).t += m_e*electron.vt + M_n*neutral.vt*dtn;
  (*momcheck).z += m_e*electron.vz + M_n*neutral.vz*dtn;

                W_0  = MAX( W_0, 1.e-14);
                IW_0 = 1./W_0;

                W2 = W2_0 - E_th;
                W  = sqrt(W2);                                 // Absolutnaya velichina otn. skorosti

                delta_v_rel.r = v_rel.r*(W*IW_0 -1.);                       // After energy loss
                delta_v_rel.t = v_rel.t*(W*IW_0 -1.);
                delta_v_rel.z = v_rel.z*(W*IW_0 -1.);

                v_rel.r += delta_v_rel.r;
                v_rel.t += delta_v_rel.t;
                v_rel.z += delta_v_rel.z;

                W  = MAX( W, 1.e-14);                         // A chtob na 0 ne delit'
                IW = 1./W;

                v_proj_rt = sqrt( SQU(v_rel.r) + SQU(v_rel.t) );     // proekciya na XY
                v_proj_rt = MAX( v_proj_rt, 1.e-14);              // ta zhe fignia
                ivp       = 1./v_proj_rt;

                cos_theta = v_rel.z * IW;               // theta - ugol mezhdu v_rel i OZ
                sin_theta = v_proj_rt*IW;

                cos_beta = v_rel.r*ivp;                // beta - ugol mezhdu proekciey v_rel na XY i os'yu abciss
                sin_beta = v_rel.t*ivp;

                 /* Sampling scattering angles in C.M. system */
                 cos_phi = 2.*RAND -1. ; // just 4 test!! Use = (2.+E-2.*pow((1+E), RAND))/E instead
                 sin_phi = sqrt(1. - SQU(cos_phi));
                 psi  = TWOPI*RAND;

                 vz = W*cos_phi;              //   komponenty posle stolknoveniya v povernutoy coor. sys.
                 vr = W*sin_phi*cos(psi);
                 vt = W*sin_phi*sin(psi);

                 v_rel.r -= vr*cos_beta*cos_theta -vt*sin_beta + vz*cos_beta*sin_theta;    // eto  uzhe minus prirascheniya otn. skor.
                 v_rel.t -= vr*sin_beta*cos_theta +vt*cos_beta + vz*sin_beta*sin_theta;    // a ne otn. skorost'
                 v_rel.z -= vz*cos_theta-vr*sin_theta;

                 v_rel.r -= delta_v_rel.r;
                 v_rel.t -= delta_v_rel.t;
                 v_rel.z -= delta_v_rel.z;

                   target_v.r = neutral.vr*dtn;
                   target_v.t = neutral.vt*dtn;
                   target_v.z = neutral.vz*dtn;

                 electron.vr -= M_n/(M_n + m_e)*v_rel.r;
                 electron.vt -= M_n/(M_n + m_e)*v_rel.t;
                 electron.vz -= M_n/(M_n + m_e)*v_rel.z;

                 target_v.r    += m_e/(M_n + m_e)*v_rel.r;   // a nado li
                 target_v.t    += m_e/(M_n + m_e)*v_rel.t;   //
                 target_v.z    += m_e/(M_n + m_e)*v_rel.z;    //               Dlia diagnostiki - prisvoy neutralam sho- nibud' nelepoe
		Particle ion;
                //Ion is born at position of neutral
                //IMPROTANT: RESULTING PARTICLE NOT SUPPOSED TO BE STORED!
                //--> in this case radical oxygen, which is not tracked in the code!
                 ion.assign_new_id();
		 ion.r = neutral.r;
                 ion.z = neutral.z;

                 ion.vr = target_v.r*dt_ion;                 // @@ scale it back
                 ion.vt = target_v.t*dt_ion;
                 ion.vz = target_v.z*dt_ion;
                 
#if DIAG_COLL
		 ion.coll_r = 0.;
		 ion.coll_t = nstep;    
#endif

*engcheck -= 0.5*(M_n-m_e)*(SQU(ion.vr*dti) + SQU(ion.vt*dti) + SQU(ion.vz*dti));

(*momcheck).r -= (M_n-m_e)*ion.vr*dti;
(*momcheck).t -= (M_n-m_e)*ion.vt*dti;
(*momcheck).z -= (M_n-m_e)*ion.vz*dti;

 /******  electron-neutral inelastic collision with a loss of E_th energy
           and birth of the ion  is done !
            now we'll do electron-electron elastic collission                ******/

                 v_rel.r= electron.vr -  target_v.r;      // snova relative velocity
                 v_rel.t= electron.vt -  target_v.t;
                 v_rel.z= electron.vz -  target_v.z;

                 W2 = SQU(v_rel.r) + SQU(v_rel.t) + SQU(v_rel.z);
                 W  = sqrt(W2);                                  // Absolutnaya velichina otn. skorosti
                 W  = MAX( W, 1.e-14);                         // A chtob na 0 ne delit'
                 IW = 1./W;

                 v_proj_rt = sqrt( SQU(v_rel.r) + SQU(v_rel.t) );     // proekciya na XY
                 v_proj_rt = MAX( v_proj_rt, 1.e-14);              // ta zhe fignia
                 ivp       = 1./v_proj_rt;

                 cos_theta = v_rel.z * IW;               // theta - ugol mezhdu v_rel i OZ
                 sin_theta = v_proj_rt*IW;

                 cos_beta = v_rel.r*ivp;                // beta - ugol mezhdu proekciey v_rel na XY i os'yu abciss
                 sin_beta = v_rel.t*ivp;

                 cos_phi = 2.*RAND -1. ;
                 sin_phi = sqrt(1 - SQU(cos_phi));
                 psi  = TWOPI*RAND;

                 vz = W*cos_phi;              //   komponenty posle stolknoveniya v povernutoy coor. sys.
                 vr = W*sin_phi*cos(psi);
                 vt = W*sin_phi*sin(psi);

                 v_rel.r -= vr*cos_beta*cos_theta -vt*sin_beta + vz*cos_beta*sin_theta;    // eto  uzhe minus prirascheniya otn. skor.
                 v_rel.t -= vr*sin_beta*cos_theta +vt*cos_beta + vz*sin_beta*sin_theta;    // a ne otn. skorost'
                 v_rel.z -= vz*cos_theta-vr*sin_theta;

                 electron.vr -= 0.5*v_rel.r;         // Note inversed sign due inversed delta be4
                 electron.vt -= 0.5*v_rel.t;         // Rastiapa!!!
                 electron.vz -= 0.5*v_rel.z;

                 target_v.r    += 0.5*v_rel.r;
                 target_v.t    += 0.5*v_rel.t;
                 target_v.z    += 0.5*v_rel.z;


 //if( *ne + inzd +1 > NE )   fprintf(stderr,"Slishkom mnogo ponastalkivali electronov");
		Particle new_electron;                	
		
		new_electron.assign_new_id();
                new_electron.r = neutral.r;
                new_electron.z = neutral.z;
                new_electron.vr = target_v.r;
                new_electron.vt = target_v.t;
                new_electron.vz = target_v.z;
#if DIAG_COLL
                 new_electron.coll_r = 0.;
		 new_electron.coll_t = nstep;
#endif
		 store_particle_properties(&new_electron, IONIZATION_ORIGIN);
		 electron_layer.store(new_electron);


 *engcheck -= 0.5*m_e*(SQU(new_electron.vr) + SQU(new_electron.vt) + SQU(new_electron.vz));
 *engcheck -= 0.5*m_e*(SQU(new_electron.vr) + SQU(new_electron.vt) + SQU(new_electron.vz));

(*momcheck).r -= m_e*(new_electron.vr + new_electron.vr);
(*momcheck).t -= m_e*(new_electron.vt + new_electron.vt);
(*momcheck).z -= m_e*(new_electron.vz + new_electron.vz);

 
	 //collision diagnostics
#if DEBUG_COLL
#if DIAG_COLL
	ncoll_ioniz[i*electron_layer.z_dim +j]++;

	diag_collisions(electron, ecoll_dr, ecoll_dt);
	diag_collisions(neutral, ncoll_dr, ncoll_dt);
#endif
#endif

	unsigned int id=k;
	neutral_cell.remove(id);
            }   //  if (W2_0 >= E_th && RAND < neutral_cell.size()*W_0*Si)
          }  //for ( i_el = 0; i_el < n_el; i_el++)
        }                    //if ( neutral_cell.size() && n_el )
      }                    //for (j=0; j < nz;  j++)
    }                     //for (i=0; i < nr;  i++)
  }
 
      
 /***********************************************************************************

  Function:     coll_el_diss_recomb_fake(...)
  Action:       electrons dissociative recombination collission (e + O2+ -> O + O )
                
Note: adapted from diss.attach, O2+ are denoted as neutrals. O atoms arent followed in this version
  
***********************************************************************************/
 
void coll_el_diss_recomb_fake( GridLayer&  neutral_layer,
                          GridLayer& electron_layer,
			  Reaction React,
                          Vec3d *momcheck, double *engcheck, int ncoll_ioniz[],
			  double ecoll_dr[], double ecoll_dt[], double ncoll_dr[], double ncoll_dt[] )
{

 const double m_e =1.;

 double  dti, dtn, VolFact;
 double  W2, W, IW, W2_0,  W_0, IW_0, v_proj_rt, ivp, cos_theta, sin_theta, cos_beta, sin_beta;
 double  cos_phi, sin_phi, psi;

 double vr, vt, vz;
 
 const double M_n=neutral_layer.mass;
 const int dt_ntrl=neutral_layer.dt_subcycling;

 Vec3d 	v_rel, delta_v_rel, target_v, v_CM;
 
 double E_th, Emin, Emax, Estep, S_i;
 int    Epoints, Eind;


  E_th  = React.Eth;
  Emin  = React.Emin;
  Emax  = React.Emax;
  Estep = React.Estep;
  Epoints = React.N;
  S_i = 0.;
  
  dti  = 1./dt_ion;      // Scaling factor 4 ion velocity
  
  dtn  = 1./dt_ntrl;

 for (int i=0; i < neutral_layer.r_dim;  ++i)
  {
   if ( i != 0 )   VolFact = 1./i;  
   else            VolFact = 6.;
     
   for (int j=0; j < neutral_layer.z_dim;  ++j)
    {
	      
	auto& neutral_cell=neutral_layer.get_cell(i,j);
	auto& electron_cell=electron_layer.get_cell(i,j);

     if ( neutral_cell.size() && electron_cell.size() )
      { 
 	   if( neutral_cell.size() > NNcell )   fprintf(stderr,"To much neutrals !!! increase NNcell");   
	    
        for (unsigned int i_el = 0; i_el < electron_cell.size() && neutral_cell.size() ; ++i_el)    // check it out
	   {
	     int k = (int)(neutral_cell.size()*RAND);
	     if (k == neutral_cell.size()) --k;              

	Particle& electron = electron_cell.particles[i_el];
	Particle& neutral = neutral_cell.particles[k];

	    v_rel.r= electron.vr - neutral.vr*dtn;      // relative velocity
	    v_rel.t= electron.vt - neutral.vt*dtn;                        // !!3.09.01 Check it on scaling 4 ion velocity
	    v_rel.z= electron.vz - neutral.vz*dtn;

	    W2_0 = SQU(v_rel.r) + SQU(v_rel.t) + SQU(v_rel.z);                // before scattering
	    W_0  = sqrt(W2_0);
   // Linear fit for Cross-section
	       if ( W2_0 >= E_th )
		{
		 Eind = (int)((W2_0 - Emin)/Estep);
		 if ( Eind < 0)
		     S_i = React.CS[0];
		 else if ( Eind >= Epoints)
		     S_i = React.CS[Epoints];
		 else
		     S_i = React.CS[Eind] + (React.CS[Eind +1] - React.CS[Eind])*((W2_0 - Emin)/Estep - Eind);
		}

		S_i *=  VolFact;

	       if (W2_0 >= E_th && RAND < neutral_cell.size()*W_0*S_i)      // !! Est' collision!! 1-exp(-neutral_cell.size()*W_0*Si*dt_coll)
		 {
   


  		   
   *engcheck += 0.5*m_e*(SQU(electron.vr) + SQU(electron.vt) + SQU(electron.vz));
   *engcheck += 0.5*M_n*(SQU(neutral.vr*dtn) + SQU(neutral.vt*dtn) + SQU(neutral.vz*dtn)) - 0.5*m_e*M_n/(m_e+M_n)*E_th;

    (*momcheck).r += m_e*electron.vr + M_n*neutral.vr*dtn;
    (*momcheck).t += m_e*electron.vt + M_n*neutral.vt*dtn;
    (*momcheck).z += m_e*electron.vz + M_n*neutral.vz*dtn;
/*
	v_CM.r += (m_e*electrons[i_el + Next_el].vr + M_n*neutrals[i_n + Next_n].vr*dtn)/(M_n + m_e);
	v_CM.t += (m_e*electrons[i_el + Next_el].vt + M_n*neutrals[i_n + Next_n].vt*dtn)/(M_n+ m_e);
	v_CM.z += (m_e*electrons[i_el + Next_el].vz + M_n*neutrals[i_n + Next_n].vz*dtn)/(M_n + m_e);


		  W2 = (W2_0 - E_th)*m_e*M_n*(M_i+M_a)/((M_n + m_e)*M_i*M_a);
		  W  = sqrt(W2);                                 // Absolutnaya velichina otn. skorosti O i O-

		  cos_phi = 2.*RAND -1. ; // just 4 test!! Use = (2.+E-2.*pow((1+E), RAND))/E instead
		  sin_phi = sqrt(1. - SQU(cos_phi));
		  psi  = TWOPI*RAND;


		  v_rel.r = W*sin_phi*cos(psi);
		  v_rel.t = W*sin_phi*sin(psi);
		  v_rel.z = W*cos_phi;

		  
		  ions[*ni + inzd].r= electrons[i_el + Next_el].r;
		  ions[*ni + inzd].z= electrons[i_el + Next_el].z;

		  ions[*ni + inzd].vr = (v_CM.r + M_a/(M_a + M_i)*v_rel.r)*dti;
		  ions[*ni + inzd].vt = (v_CM.t + M_a/(M_a + M_i)*v_rel.t)*dti;
		  ions[*ni + inzd].vz = (v_CM.z + M_a/(M_a + M_i)*v_rel.z)*dti;

#if DIAG_COLL
		  ions[*ni + inzd].coll_r = 0.;
		  ions[*ni + inzd].coll_t = nstep;    
#endif

		  
		  atoms[*na + inzd].r= electrons[i_el + Next_el].r;
		  atoms[*na + inzd].z= electrons[i_el + Next_el].z;

		  atoms[*na + inzd].vr = (v_CM.r - M_i/(M_a + M_i)*v_rel.r)*dti;
		  atoms[*na + inzd].vt = (v_CM.t - M_i/(M_a + M_i)*v_rel.t)*dti;
		  atoms[*na + inzd].vz = (v_CM.z - M_i/(M_a + M_i)*v_rel.z)*dti;

#if DIAG_COLL
		  atoms[*na + inzd].coll_r = 0.;
		  atoms[*na + inzd].coll_t = nstep;
#endif



  *engcheck -= 0.5*(M_i)*(SQU(ions[*ni + inzd].vr*dti) + SQU(ions[*ni + inzd].vt*dti) + SQU(ions[*ni + inzd].vz*dti));
  *engcheck -= 0.5*(M_a)*(SQU(atoms[*na + inzd].vr*dti) + SQU(atoms[*na + inzd].vt*dti) + SQU(atoms[*na + inzd].vz*dti));

  (*momcheck).r -= (M_i*ions[*ni + inzd].vr + M_a*atoms[*na + inzd].vr)*dti;
  (*momcheck).t -= (M_i*ions[*ni + inzd].vt + M_a*atoms[*na + inzd].vt)*dti;
  (*momcheck).z -= (M_i*ions[*ni + inzd].vz + M_a*atoms[*na + inzd].vz)*dti;

		    ++inzd;
		    --N2coll;
   */
	   //collision diagnostics
#if DEBUG_COLL
#if DIAG_COLL
	  ncoll_ioniz[i*electron_layer.z_dim + j]++;
	diag_collisions(electron, ecoll_dr, ecoll_dt);
	diag_collisions(neutral, ncoll_dr, ncoll_dt);
	  
#endif
#endif
	unsigned int id=k;
	neutral_cell.remove(id);
	electron_cell.remove(i_el);
	  
	      }   //  if (W2_0 >= E_th && RAND < neutral_cell.size()*W_0*Si)
	    }  //for ( i_el = 0; i_el < n_el; i_el++)
          }//if(n_el<neutral_cell.size())
     }                    //for (j=0; j < nz;  j++)
    }                     //for (i=0; i < nr;  i++)
  }
 
 /***********************************************************************************

  Function:     coll_ip_im_neuutralization_fake(...)
  Action:       positive -negative ion dissociative collision (O2+  + O- -> O2 + O )
  
  neutrals[] - positive ions
  ions[] - negative ions

Note: was derived from diss attach, produced neutrals arent followed
***********************************************************************************/
 
void coll_ip_im_neutralization_fake( GridLayer&  neutral_layer,
                          GridLayer& ion_layer,
			  Reaction React,
                          Vec3d *momcheck, double *engcheck, int ncoll_ioniz[],
			  double ecoll_dr[], double ecoll_dt[], double ncoll_dr[], double ncoll_dt[] )
{

 const double m_e =1.;

 double  dti, dtn, VolFact;
 double  W2, W, IW, W2_0,  W_0, IW_0, v_proj_rt, ivp, cos_theta, sin_theta, cos_beta, sin_beta;
 double  cos_phi, sin_phi, psi;

 double vr, vt, vz;
 
 const double M_n=neutral_layer.mass;
 const double M_i=ion_layer.mass;
 const int dt_ntrl=neutral_layer.dt_subcycling;
 const int dt_ion=ion_layer.dt_subcycling;

 Vec3d 	v_rel, delta_v_rel, target_v, v_CM;
 
 double E_th, Emin, Emax, Estep, S_i;
 int    Epoints, Eind;


  E_th  = React.Eth;
  Emin  = React.Emin;
  Emax  = React.Emax;
  Estep = React.Estep;
  Epoints = React.N;
  S_i = 0.;
  

  dti  = 1./dt_ion;      // Scaling factor 4 ion velocity
  
  dtn  = 1./dt_ntrl;

 for (int i=0; i < neutral_layer.r_dim;  ++i)
  {
   if ( i != 0 )   VolFact = 1./i;  
   else            VolFact = 6.;
     
   for (int j=0; j < neutral_layer.z_dim;  ++j)
    {

	auto& neutral_cell=neutral_layer.get_cell(i,j);
	auto& ion_cell=ion_layer.get_cell(i,j);


     if ( neutral_cell.size() && ion_cell.size() )
      { 
	    if( neutral_cell.size() > NNcell )   fprintf(stderr,"To much neutrals !!! increase NNcell");   
	    
        for (unsigned int i_i = 0; i_i < ion_cell.size() && neutral_cell.size() ; ++i_i)    // check it out
         {
           int k = (int)(neutral_cell.size()*RAND);
           if (k == neutral_cell.size()) --k;              

	Particle neutral = neutral_cell.particles[k];		
	Particle ion = ion_cell.particles[i_i];		

          v_rel.r= ion.vr - neutral.vr*dtn;      // relative velocity
          v_rel.t= ion.vt - neutral.vt*dtn;                        // !!3.09.01 Check it on scaling 4 ion velocity
          v_rel.z= ion.vz - neutral.vz*dtn;

          W2_0 = SQU(v_rel.r) + SQU(v_rel.t) + SQU(v_rel.z);                // before scattering
          W_0  = sqrt(W2_0);
 // Linear fit for Cross-section
             if ( W2_0 >= E_th )
              {
               Eind = (int)((W2_0 - Emin)/Estep);
               if ( Eind < 0)
                   S_i = React.CS[0];
               else if ( Eind >= Epoints)
                   S_i = React.CS[Epoints];
               else
                   S_i = React.CS[Eind] + (React.CS[Eind +1] - React.CS[Eind])*((W2_0 - Emin)/Estep - Eind);
              }

              S_i *=  VolFact;

             if (W2_0 >= E_th && RAND < neutral_cell.size()*W_0*S_i)      // !! Est' collision!! 1-exp(-neutral_cell.size()*W_0*Si*dt_coll)
               {
 
 *engcheck += 0.5*M_i*(SQU(ion.vr) + SQU(ion.vt) + SQU(ion.vz));
 *engcheck += 0.5*M_n*(SQU(neutral.vr*dtn) + SQU(neutral.vt*dtn) + SQU(neutral.vz*dtn)) - 0.5*M_i*M_n/(M_i+M_n)*E_th;

  (*momcheck).r += M_i*ion.vr + M_n*neutral.vr*dtn;
  (*momcheck).t += M_i*ion.vt + M_n*neutral.vt*dtn;
  (*momcheck).z += M_i*ion.vz + M_n*neutral.vz*dtn;

/*		   
      v_CM.r += (m_e*electrons[i_el + Next_el].vr + M_n*neutrals[i_n + Next_n].vr*dtn)/(M_n + m_e);
      v_CM.t += (m_e*electrons[i_el + Next_el].vt + M_n*neutrals[i_n + Next_n].vt*dtn)/(M_n+ m_e);
      v_CM.z += (m_e*electrons[i_el + Next_el].vz + M_n*neutrals[i_n + Next_n].vz*dtn)/(M_n + m_e);


                W2 = (W2_0 - E_th)*m_e*M_n*(M_i+M_a)/((M_n + m_e)*M_i*M_a);
                W  = sqrt(W2);                                 // Absolutnaya velichina otn. skorosti O i O-

	     	cos_phi = 2.*RAND -1. ; // just 4 test!! Use = (2.+E-2.*pow((1+E), RAND))/E instead
                sin_phi = sqrt(1. - SQU(cos_phi));
                psi  = TWOPI*RAND;


		v_rel.r = W*sin_phi*cos(psi);
                v_rel.t = W*sin_phi*sin(psi);
                v_rel.z = W*cos_phi;

		
		ions[*ni + inzd].r= electrons[i_el + Next_el].r;
		ions[*ni + inzd].z= electrons[i_el + Next_el].z;

		ions[*ni + inzd].vr = (v_CM.r + M_a/(M_a + M_i)*v_rel.r)*dti;
		ions[*ni + inzd].vt = (v_CM.t + M_a/(M_a + M_i)*v_rel.t)*dti;
		ions[*ni + inzd].vz = (v_CM.z + M_a/(M_a + M_i)*v_rel.z)*dti;

#if DIAG_COLL
		ions[*ni + inzd].coll_r = 0.;
		ions[*ni + inzd].coll_t = nstep;    
#endif

		
		atoms[*na + inzd].r= electrons[i_el + Next_el].r;
		atoms[*na + inzd].z= electrons[i_el + Next_el].z;

		atoms[*na + inzd].vr = (v_CM.r - M_i/(M_a + M_i)*v_rel.r)*dti;
		atoms[*na + inzd].vt = (v_CM.t - M_i/(M_a + M_i)*v_rel.t)*dti;
		atoms[*na + inzd].vz = (v_CM.z - M_i/(M_a + M_i)*v_rel.z)*dti;

#if DIAG_COLL
		atoms[*na + inzd].coll_r = 0.;
		atoms[*na + inzd].coll_t = nstep;
#endif



*engcheck -= 0.5*(M_i)*(SQU(ions[*ni + inzd].vr*dti) + SQU(ions[*ni + inzd].vt*dti) + SQU(ions[*ni + inzd].vz*dti));
*engcheck -= 0.5*(M_a)*(SQU(atoms[*na + inzd].vr*dti) + SQU(atoms[*na + inzd].vt*dti) + SQU(atoms[*na + inzd].vz*dti));

(*momcheck).r -= (M_i*ions[*ni + inzd].vr + M_a*atoms[*na + inzd].vr)*dti;
(*momcheck).t -= (M_i*ions[*ni + inzd].vt + M_a*atoms[*na + inzd].vt)*dti;
(*momcheck).z -= (M_i*ions[*ni + inzd].vz + M_a*atoms[*na + inzd].vz)*dti;

                  ++inzd;
                  --N2coll;
 */
	 //collision diagnostics
#if DEBUG_COLL
#if DIAG_COLL
	ncoll_ioniz[i*neutral_layer.z_dim + j]++;
	diag_collisions(neutral, ncoll_dr, ncoll_dt);
	diag_collisions(ion, icoll_dr, icoll_dt);


#endif
#endif

	unsigned int id = k;
	neutral_cell.remove( id );
	ion_cell.remove( i_i );

//TODO add created O2 if !NTRL_CONST

		
            }   //  if (W2_0 >= E_th && RAND < neutral_cell.size()*W_0*Si)
          }  //for ( i_el = 0; i_el < n_el; i_el++)
	}
      }                    //for (j=0; j < nz;  j++)
    }                     //for (i=0; i < nr;  i++)
  }
  
/***********************************************************************************

  Function:     coll_neutr_detachment_fake(...)
  Action:       neutrals detachment  (O- + O2 -> O + O2 +e )
              
  neutrals are O2
  ions are O-
  we dont follow produced neutrals, only charged particles are followed
  
***********************************************************************************/
 

 void coll_neutr_detachment_fake( GridLayer&  neutral_layer,
                          GridLayer& ion_layer,
                          GridLayer& electron_layer,
			  Reaction React,
                          Vec3d *momcheck, double *engcheck, int ncoll_ioniz[],
			  double ecoll_dr[], double ecoll_dt[], double ncoll_dr[], double ncoll_dt[] )
{

 const double m_e =1.;

 double  dti, dtn, dti_over_dtn, VolFact;
 double  W2, W,  W2_0,  W_0;
 double  cos_phi, sin_phi, psi;

 double vr, vt, vz;
 
 const double M_n=neutral_layer.mass;
 const double M_i=ion_layer.mass;
 const int dt_ntrl=neutral_layer.dt_subcycling;
 const int dt_ion=ion_layer.dt_subcycling;

 Vec3d 	v_rel, delta_v_rel, target_v;
 
 double E_th, Emin, Emax, Estep, S_i;
 int    Epoints, Eind;


  E_th  = React.Eth;
  Emin  = React.Emin;
  Emax  = React.Emax;
  Estep = React.Estep;
  Epoints = React.N;
  S_i = 0.;
  
  dti  = 1./dt_ion;      // Scaling factor 4 ion velocity
  
  dtn  = 1./dt_ntrl;

  dti_over_dtn=dti/dtn;

 for (int i=0; i < neutral_layer.r_dim;  ++i)
  {
   if ( i != 0 )   VolFact = 1./i;  
   else            VolFact = 6.;
     
   for (int j=0; j < neutral_layer.z_dim;  ++j)
    {
	auto& neutral_cell=neutral_layer.get_cell(i,j);
	auto& ion_cell=ion_layer.get_cell(i,j);
	auto& electron_cell=electron_layer.get_cell(i,j);

     if ( neutral_cell.size() && ion_layer.size() )
      { 
	    if( neutral_cell.size() > NNcell )   fprintf(stderr,"To much neutrals !!! increase NNcell");   
	    
        for (unsigned int i_i = 0; i_i < ion_cell.size() && neutral_cell.size() ; ++i_i)    // check it out
         {
           int k = (int)(neutral_cell.size()*RAND);
           if (k == neutral_cell.size()) --k;              

	Particle neutral = neutral_cell.particles[k];		
	Particle ion = ion_cell.particles[i_i];		

          v_rel.r= ion.vr*dti_over_dtn - neutral.vr;      // relative velocity
          v_rel.t= ion.vt*dti_over_dtn - neutral.vt;                        // !!3.09.01 Check it on scaling 4 ion velocity
          v_rel.z= ion.vz*dti_over_dtn - neutral.vz;

          W2_0 = SQU(v_rel.r) + SQU(v_rel.t) + SQU(v_rel.z);                // before scattering
          W_0  = sqrt(W2_0);
 // Linear fit for Cross-section
             if ( W2_0 >= E_th )
              {
               Eind = (int)((W2_0 - Emin)/Estep);
               if ( Eind < 0)
                   S_i = React.CS[0];
               else if ( Eind >= Epoints)
                   S_i = React.CS[Epoints];
               else
                   S_i = React.CS[Eind] + (React.CS[Eind +1] - React.CS[Eind])*((W2_0 - Emin)/Estep - Eind);
              }

              S_i *=  VolFact;

             if (W2_0 >= E_th && RAND < neutral_cell.size()*W_0*S_i)      // !! Est' collision!! 1-exp(-neutral_cell.size()*W_0*Si*dt_coll)
               {
 
		   /*
 *engcheck += 0.5*m_e*(SQU(electrons[i_el + Next_el].vr) + SQU(electrons[i_el + Next_el].vt) + SQU(electrons[i_el + Next_el].vz));
 *engcheck += 0.5*M_n*(SQU(neutrals[i_n + Next_n].vr*dtn) + SQU(neutrals[i_n + Next_n].vt*dtn) + SQU(neutrals[i_n + Next_n].vz*dtn)) - 0.5*m_e*M_n/(m_e+M_n)*E_th;

  (*momcheck).r += m_e*electrons[i_el + Next_el].vr + M_n*neutrals[i_n + Next_n].vr*dtn;
  (*momcheck).t += m_e*electrons[i_el + Next_el].vt + M_n*neutrals[i_n + Next_n].vt*dtn;
  (*momcheck).z += m_e*electrons[i_el + Next_el].vz + M_n*neutrals[i_n + Next_n].vz*dtn;
*/
                W2 = 0.1*( W2_0 - E_th)*M_i*M_n/((M_n + M_i)*m_e);
                W  = sqrt(W2);                                 // Absolutnaya velichina otn. skorosti

                           /* Sampling scattering angles in C.M. system */
                 cos_phi = 2.*RAND -1. ; // just 4 test!! Use = (2.+E-2.*pow((1+E), RAND))/E instead
                 sin_phi = sqrt(1. - SQU(cos_phi));
                 psi  = TWOPI*RAND;

                 v_rel.r = W*cos_phi;              //   komponenty posle stolknoveniya v povernutoy coor. sys.
                 v_rel.t = W*sin_phi*cos(psi);
                 v_rel.z = W*sin_phi*sin(psi);

		//electron is born
		Particle electron;
		electron.assign_new_id();
		
                 electron.r = ion.r;
                 electron.z = ion.z;

		 electron.vr = v_rel.r;
                 electron.vt = v_rel.t;
                 electron.vz = v_rel.z;
 
#if DIAG_COLL
		 electron.coll_r = 0.;
		 electron.coll_t = nstep;    
#endif
 		store_particle_properties (&electron, IONIZATION_ORIGIN);

		electron_layer.store(electron);

	 //collision diagnostics
	 //
#if DEBUG_COLL
#if DIAG_COLL
	ncoll_ioniz[i*neutral_layer.z_dim + j]++;
	diag_collisions(neutral, ncoll_dr, ncoll_dt);
	diag_collisions(ion, icoll_dr, icoll_dt);
#endif
#endif
	ion_cell.remove( i_i ); 	    

           }   //  if (W2_0 >= E_th && RAND < neutral_cell.size()*W_0*Si)
         }  //for ( i_el = 0; i_el < n_el; i_el++)
       }                    //if ( neutral_cell.size() && n_el )
      }                    //for (j=0; j < nz;  j++)
    }                     //for (i=0; i < nr;  i++)
  }


