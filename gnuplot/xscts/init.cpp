void init_reactions_Xe(void) {
  double E;
  double tmp;
  double Ecoeff, CScoeff;
  double Emin, Emax, Eth;

  FILE *outputfile;

  double edata[80], Sdata[80];

  int Nmax;

  Nmax = Npoints_CS - 1;


  /**************************************************************************************/
  /*                            Xe + e -> Xe+ + e 					*/
  /*                         Ionization cross-section 					*/
  /*                    Source:  Hayashi NIFS report   2003 				*/
  /*                              http://www.nifs.ac.jp/report/NIFS-DATA-079.pdf p.134  */
  /**************************************************************************************/
  //    E <eV>            S <10^-16 cm2>
  edata[0] = 12.13;  	Sdata[0] = 0.;
  edata[1] = 12.5;  	Sdata[1] = 0.105;
  edata[2] = 13.;  	Sdata[2] = 0.247;
  edata[3] = 13.5;  	Sdata[3] = 0.409;
  edata[4] = 14.;  	Sdata[4] = 0.57;
  edata[5] = 14.5;  	Sdata[5] = 0.751;
  edata[6] = 15.;  	Sdata[6] = 0.931;
  edata[7] = 15.5;  	Sdata[7] = 1.09;
  edata[8] = 16.;  	Sdata[8] = 1.25;
  edata[9] = 17.;  	Sdata[9] = 1.56;
  edata[10] = 18.;  	Sdata[10] = 1.82;
  edata[11] = 19.;  	Sdata[11] = 2.06;
  edata[12] = 20.;  	Sdata[12] = 2.30;
  edata[13] = 22.;  	Sdata[13] = 2.76;
  edata[14] = 24.;  	Sdata[14] = 3.14;
  edata[15] = 28.;  	Sdata[15] = 3.70;
  edata[16] = 31.;  	Sdata[16] = 4.02;
  edata[17] = 36.;  	Sdata[17] = 4.43;
  edata[18] = 45.;  	Sdata[18] = 4.87;
  edata[19] = 60.;  	Sdata[19] = 5.34;
  edata[20] = 100.;  	Sdata[20] = 5.76;
  edata[21] = 110.;  	Sdata[21] = 5.88;
  edata[22] = 125.;  	Sdata[22] = 5.88;
  edata[23] = 135.;  	Sdata[23] = 5.8;
  edata[24] = 145.;  	Sdata[24] = 5.71;
  edata[25] = 170.;  	Sdata[25] = 5.36;
  edata[26] = 200.;  	Sdata[26] = 5.02;
  edata[27] = 250.;  	Sdata[27] = 4.37;
  edata[28] = 300.;  	Sdata[28] = 3.96;
  edata[29] = 500.;  	Sdata[29] = 2.97;
  edata[30] = 700.;  	Sdata[30] = 2.43;
  edata[31] = 1000.;  	Sdata[31] = 1.90;
  edata[32] = 10000.;  	Sdata[32] = 0.0;

  Emin = 12.13;  // eV
  Eth = 12.13;   // eV
  Emax = 200.;   // eV

  Ecoeff = (dt / dr) * (dt / dr) * (1. + mn_over_me) / mn_over_me / (0.5 * T_e0);

  CScoeff = collision_amplification_factor * ScaleF * n_e0 * dr_0 * ncoll_n / Ncell1;  // Amplification factor was  50.

  React_Xe_i.N = Nmax;
  React_Xe_i.Emin = Emin * Ecoeff;
  React_Xe_i.Eth = Eth * Ecoeff;
  React_Xe_i.Emax = Emax * Ecoeff;
  React_Xe_i.Estep = (React_Xe_i.Emax - React_Xe_i.Emin) / React_Xe_i.N;

  outputfile = fopen("Xe_i.dat", "w");

  for (int i = 0; i <= React_Xe_i.N; i++) {
    E = i * React_Xe_i.Estep + React_Xe_i.Emin;

    tmp = 0.;

    if (E / Ecoeff <= edata[0])
      tmp = Sdata[0];
    else if (E / Ecoeff >= edata[26])
      tmp = Sdata[26];
    else {
      int ind = 0;
      while (E / Ecoeff > edata[ind]) ind++;
      assert( ind >= 1 && "calculated index cases out of bounds access to Sdata" );
      tmp = Sdata[ind - 1] +
	    (Sdata[ind] - Sdata[ind - 1]) * (E / Ecoeff - edata[ind - 1]) /
		(edata[ind] - edata[ind - 1]);
    }

    React_Xe_i.CS[i] = CScoeff * tmp * 1.e-16 * coll_fac_ntrl_ntrl;   // [1]

    fprintf(outputfile, "%f  %15.7e \n", E / Ecoeff, React_Xe_i.CS[i] / CScoeff); // [eV] [cm^2] (scaled)
  }
  fclose(outputfile);


  /**************************************************************************************/
  /*                             Xe + e -> Xe + e 					*/
  /*                         Total elastic cross-section 				*/
  /*                     Source:  Hayashi NIFS report   2003 			   	*/
  /*                              http://www.nifs.ac.jp/report/NIFS-DATA-079.pdf p.129  */
  /**************************************************************************************/
  //    E <eV>            S <10^-16 cm2>
  edata[0] = 0.001;  	Sdata[0] = 123.;
  edata[1] = 0.005;  	Sdata[1] = 99.0;
  edata[2] = 0.01;  	Sdata[2] = 84.3;
  edata[3] = 0.04;  	Sdata[3] = 47.2;
  edata[4] = 0.1;  	Sdata[4] = 21.5;
  edata[5] = 0.16;  	Sdata[5] = 11.4;
  edata[6] = 0.2;  	Sdata[6] = 7.79;
  edata[7] = 0.32;  	Sdata[7] = 2.87;
  edata[8] = 0.35;  	Sdata[8] = 2.16;
  edata[9] = 0.4;  	Sdata[9] = 1.47;
  edata[10] = 0.5;  	Sdata[10] = 0.64;
  edata[11] = 0.55;  	Sdata[11] = 0.45;
  edata[12] = 0.6;  	Sdata[12] = 0.34;
  edata[13] = 0.64;  	Sdata[13] = 0.315;
  edata[14] = 0.7;  	Sdata[14] = 0.36;
  edata[15] = 0.8;  	Sdata[15] = 0.54;
  edata[16] = 0.9;  	Sdata[16] = 0.81;
  edata[17] = 1.0;  	Sdata[17] = 1.21;
  edata[18] = 1.2;  	Sdata[18] = 2.21;
  edata[19] = 1.4;  	Sdata[19] = 3.38;
  edata[20] = 1.6;  	Sdata[20] = 4.68;
  edata[21] = 2.;  	Sdata[21] = 7.39;
  edata[22] = 2.6;  	Sdata[22] = 12.4;
  edata[23] = 3.0;  	Sdata[23] = 16.1;
  edata[24] = 3.6;  	Sdata[24] = 21.1;
  edata[25] = 4.4;  	Sdata[25] = 26.2;
  edata[26] = 5.0;  	Sdata[26] = 28.3;
  edata[27] = 5.4;  	Sdata[27] = 28.8;
  edata[28] = 6.0;  	Sdata[28] = 28.0;
  edata[29] = 8.0;  	Sdata[29] = 22.5;
  edata[30] = 10.;  	Sdata[30] = 17.;
  edata[31] = 12.;	Sdata[31] = 13.1;
  edata[32] = 15.;  	Sdata[32] = 9.04;
  edata[33] = 20.;  	Sdata[33] = 6.12;
  edata[34] = 30.;  	Sdata[34] = 5.11;
  edata[35] = 40.;  	Sdata[35] = 4.05;
  edata[36] = 55.;  	Sdata[36] = 2.30;
  edata[37] = 70.;  	Sdata[37] = 1.52;
  edata[38] = 100.;  	Sdata[38] = 1.50;
  edata[39] = 125.;  	Sdata[39] = 1.73;
  edata[40] = 140.;  	Sdata[40] = 1.71;
  edata[41] = 200.;  	Sdata[41] = 1.24;
  edata[42] = 400.;  	Sdata[42] = 0.8;
  edata[43] = 600.;  	Sdata[43] = 0.71;
  edata[44] = 1000.;  	Sdata[44] = 0.49;

  Emin = 0.001;  // eV
  Eth = 0.0;     // eV
  Emax = 100.;   // eV

  Ecoeff = (dt / dr) * (dt / dr) * (1. + mn_over_me) / mn_over_me / (0.5 * T_e0);

  CScoeff = collision_amplification_factor * ScaleF * n_e0 * dr_0 * ncoll_n / Ncell1;  // Amplification factor was  50.

  React_Xe_el.N = Nmax;
  React_Xe_el.Emin = Emin * Ecoeff;
  React_Xe_el.Eth = Eth * Ecoeff;
  React_Xe_el.Emax = Emax * Ecoeff;
  React_Xe_el.Estep = (React_Xe_el.Emax - React_Xe_el.Emin) / React_Xe_el.N;

  outputfile = fopen("Xe_el.dat", "w");


  for (int i = 0; i <= React_Xe_el.N; i++) {
    E = i * React_Xe_el.Estep + React_Xe_el.Emin;

    tmp = 0.;

   if (E / Ecoeff <= edata[0]) 
      tmp = Sdata[0];
    else if (E / Ecoeff >= edata[38])
      tmp = Sdata[38];
    else {
      int ind = 0;
      while (E / Ecoeff > edata[ind]) ind++;
      assert( ind >= 1 && "calculated index cases out of bounds access to Sdata" );
      tmp = Sdata[ind - 1] +
	    (Sdata[ind] - Sdata[ind - 1]) * (E / Ecoeff - edata[ind - 1]) /
		(edata[ind] - edata[ind - 1]);
    }

    React_Xe_el.CS[i] = CScoeff * tmp * 1.e-16 * coll_fac_ntrl_ntrl;  // [1]

    fprintf(outputfile, "%f  %15.7e \n", E / Ecoeff, React_Xe_el.CS[i] / CScoeff); // [eV] [cm^2] (scaled)
  }
  fclose(outputfile);


  /**************************************************************************************/
  /*                      Xe + e -> Xe* + e 						*/
  /*                Total excitation cross-section 					*/
  /*       Source:  M. Hayashi, Journal of Physics D: 16 (1983) 581 			*/
  /**************************************************************************************/
  //    E <eV>            S <10^-16 cm^2>
  edata[0] = 8.32;  	Sdata[0] = 0.0;
  edata[1] = 8.5;  	Sdata[1] = 0.026;
  edata[2] = 9.;  	Sdata[2] = 0.126;
  edata[3] = 9.5;  	Sdata[3] = 0.131;
  edata[4] = 10.;  	Sdata[4] = 0.18;
  edata[5] = 10.5;  	Sdata[5] = 0.24;
  edata[6] = 11.;  	Sdata[6] = 0.42;
  edata[7] = 11.5;  	Sdata[7] = 0.62;
  edata[8] = 12.;  	Sdata[8] = 0.84;
  edata[9] = 12.5;  	Sdata[9] = 1.05;
  edata[10] = 13.;  	Sdata[10] = 1.28;
  edata[11] = 14.;  	Sdata[11] = 1.70;
  edata[12] = 15.;  	Sdata[12] = 2.14;
  edata[13] = 16.;  	Sdata[13] = 2.55;
  edata[14] = 18.;  	Sdata[14] = 3.35;
  edata[15] = 20.;  	Sdata[15] = 3.73;
  edata[16] = 25.;  	Sdata[16] = 3.85;
  edata[17] = 30.;  	Sdata[17] = 3.57;
  edata[18] = 40.;  	Sdata[18] = 2.85;
  edata[19] = 50.;  	Sdata[19] = 2.40;
  edata[20] = 60.;  	Sdata[20] = 2.1;
  edata[21] = 70.;  	Sdata[21] = 1.85;
  edata[22] = 80.;  	Sdata[22] = 1.66;
  edata[23] = 90.;  	Sdata[23] = 1.52;
  edata[24] = 100.;  	Sdata[24] = 1.38;
  edata[25] = 150.;  	Sdata[25] = 1.0;
  edata[26] = 200.;  	Sdata[26] = 0.8;
  edata[27] = 300.;  	Sdata[27] = 0.568;
  edata[28] = 400.;  	Sdata[28] = 0.465;
  edata[29] = 500.;  	Sdata[29] = 0.395;
  edata[30] = 600.;  	Sdata[30] = 0.344;
  edata[31] = 700.;  	Sdata[31] = 0.302;
  edata[32] = 800.;  	Sdata[32] = 0.277;
  edata[33] = 900.;  	Sdata[33] = 0.252;
  edata[34] = 1000.;  	Sdata[34] = 0.231;
  edata[35] = 2000.;  	Sdata[35] = 0.132;
  edata[36] = 3000.;  	Sdata[36] = 0.095;
  edata[37] = 4000.;  	Sdata[37] = 0.075;
  edata[38] = 5000.;  	Sdata[38] = 0.063;
  edata[39] = 10000.;  	Sdata[39] = 0.036;

  Emin = 8.32;  // eV
  Eth = 8.32;  // eV
  Emax = 500.;  // eV

  Ecoeff = (dt / dr) * (dt / dr) * (1. + mn_over_me) / mn_over_me / (0.5 * T_e0);

  CScoeff = collision_amplification_factor * ScaleF * n_e0 * dr_0 * ncoll_n / Ncell1;  
  
  React_Xe_tex.N = Nmax;
  React_Xe_tex.Emin = Emin * Ecoeff;
  React_Xe_tex.Eth = Eth * Ecoeff;
  React_Xe_tex.Emax = Emax * Ecoeff;
  React_Xe_tex.Estep = (React_Xe_tex.Emax - React_Xe_tex.Emin) / React_Xe_tex.N;

  outputfile = fopen("Xe_tex.dat", "w");

  for (int i = 0; i <= React_Xe_tex.N; i++) {
    E = i * React_Xe_tex.Estep + React_Xe_tex.Emin;

    tmp = 0.;

    if (E / Ecoeff <= edata[0])
      tmp = Sdata[0];
    else if (E / Ecoeff >= edata[29])
      tmp = Sdata[29];
    else {
      int ind = 0;
      while (E / Ecoeff > edata[ind]) ind++;
      assert( ind >= 1 && "calculated index cases out of bounds access to Sdata" );
      tmp = Sdata[ind - 1] +
	    (Sdata[ind] - Sdata[ind - 1]) * (E / Ecoeff - edata[ind - 1]) /
		(edata[ind] - edata[ind - 1]);
    }

    React_Xe_tex.CS[i] = CScoeff * tmp * 1.e-16 * coll_fac_ntrl_ntrl;  // [1]

    fprintf(outputfile, "%f  %15.7e \n", E / Ecoeff, React_Xe_tex.CS[i] / CScoeff); // [eV] [cm^2] (scaled)
  }
  fclose(outputfile);


  /********************************************************************************/
  /*                      Xe + Xe+ -> Xe + Xe+ 					  */
  /*        Elastic and Charge Exchange  (isotropic scattering aproximation) 	  */
  /* Source: http://jila.colorado.edu/~avp/collision_data/ionneutral/IONATOM.TXT  */
  /********************************************************************************/
  Emin = 0.01;    // eV
  Eth = 0.0;      // eV
  Emax = 300.01;  // eV

  Ecoeff = (dt / dr) * (dt / dr) * (mi_over_me + mn_over_me) / mn_over_me /
	   mi_over_me / (0.5 * T_e0) * dt_ntrl * dt_ntrl; 

  CScoeff = collision_amplification_factor * ScaleF * n_e0 * dr_0 * ncoll_n / Ncell1 / dt_ntrl;  

  React_Xep_Xe.N = Nmax;
  React_Xep_Xe.Emin = Emin * Ecoeff;
  React_Xep_Xe.Eth = Eth * Ecoeff;
  React_Xep_Xe.Emax = Emax * Ecoeff;
  React_Xep_Xe.Estep = (React_Xep_Xe.Emax - React_Xep_Xe.Emin) / React_Xep_Xe.N;

  outputfile = fopen("Xep_Xe.dat", "w");

  for (int i = 0; i <= React_Xep_Xe.N; i++) {
    E = i * React_Xep_Xe.Estep + React_Xep_Xe.Emin;

    double Eev = E / Ecoeff;

    tmp = 3.39e-15 / sqrt(Eev) +
	  7.2e-15 / pow(Eev, 0.42) * pow(1. + pow(10. * Eev, 2.), 0.2) /
	      (1. + pow(0.09 / Eev, 1.3)) / (1. + pow(0.001 * Eev, 0.25));

    if(!old_Xenon_CX_CS){
      React_Xep_Xe.CS[i] = CScoeff * tmp * coll_fac_ntrl_ntrl;   // correct collisions!!! need new HEMP runs :/
    }else{
      React_Xep_Xe.CS[i] = CScoeff * tmp / dt_ion * coll_fac_ntrl_ntrl;  // old CS for continuing old runs!! :/
      if(i==0){ print0("!!!! Wrongly scaled CX-collisions cross section for Xenon !!!! (as used in old runs)\n"); }
    }

    fprintf(outputfile, "%f  %15.7e \n", E / Ecoeff, React_Xep_Xe.CS[i] / CScoeff); // [eV] [cm^2] (scaled)
  }
  fclose(outputfile);


  /***************************************************************************************/
  /*                      Xe + Xe ->  Xe + Xe 						 */
  /*                Atom-Atom elastic collision 					 */
  /*   Source:  http://jila.colorado.edu/~avp/collision_data/neutralneutral/atomatom.txt */
  /***************************************************************************************/
  // edata in [eV]      and     Sdata in [m^2]
  edata[0] = 1.00000e-02;  Sdata[0] = 9.96731e-18;
  edata[1] = 1.77828e-02;  Sdata[1] = 8.59084e-18;
  edata[2] = 3.16228e-02;  Sdata[2] = 7.84381e-18;
  edata[3] = 5.62341e-02;  Sdata[3] = 7.17105e-18;
  edata[4] = 1.00000e-01;  Sdata[4] = 6.35342e-18;
  edata[5] = 1.77828e-01;  Sdata[5] = 5.51434e-18;
  edata[6] = 3.16228e-01;  Sdata[6] = 4.88312e-18;
  edata[7] = 5.62341e-01;  Sdata[7] = 4.33354e-18;
  edata[8] = 1.00000e+00;  Sdata[8] = 3.98392e-18;
  edata[9] = 1.77828e+00;  Sdata[9] = 3.79373e-18;
  edata[10] = 3.16228e+00;  Sdata[10] = 3.17573e-18;
  edata[11] = 5.62341e+00;  Sdata[11] = 2.79482e-18;
  edata[12] = 1.00000e+01;  Sdata[12] = 2.92242e-18;
  edata[13] = 1.77828e+01;  Sdata[13] = 2.70825e-18;
  edata[14] = 3.16228e+01;  Sdata[14] = 2.22229e-18;
  edata[15] = 5.62341e+01;  Sdata[15] = 1.74008e-18;
  edata[16] = 1.00000e+02;  Sdata[16] = 1.37795e-18;
  edata[17] = 1.77828e+02;  Sdata[17] = 1.13733e-18;
  edata[18] = 3.16228e+02;  Sdata[18] = 9.85212e-19;
  edata[19] = 5.62341e+02;  Sdata[19] = 8.89419e-19;
  edata[20] = 1.00000e+03;  Sdata[20] = 8.26278e-19;
  edata[21] = 1.77828e+03;  Sdata[21] = 7.81108e-19;
  edata[22] = 3.16228e+03;  Sdata[22] = 7.45181e-19;
  edata[23] = 5.62341e+03;  Sdata[23] = 7.13897e-19;

  Emin = 1.0e-2;      // eV
  // Eth =  0.0;      // eV
  Emax = 3.16228e+02; // eV

  Ecoeff = (dt / dr) * (dt / dr) * 2. / mn_over_me / (0.5 * T_e0) * dt_ntrl * dt_ntrl;
  
  CScoeff = collision_amplification_factor * ScaleF * n_e0 * dr_0 * ncoll_n_n / Ncell1 / dt_ntrl; 


  React_Xe_Xe.N = Nmax;
  React_Xe_Xe.Emin = Emin * Ecoeff;

  React_Xe_Xe.Emax = Emax * Ecoeff;
  React_Xe_Xe.Estep = (React_Xe_Xe.Emax - React_Xe_Xe.Emin) / React_Xe_Xe.N;
  React_Xe_Xe.Eth = sqrt(edata[18] * Ecoeff) * Sdata[18] * 1.e4 * CScoeff; // formular from which source??
  //printf("-> Xe-Xe MAX(Sigma*v) = %f \n",  React_Xe_Xe.Eth );

  outputfile = fopen("Xe_Xe.dat", "w");

  for (int i = 0; i <= React_Xe_Xe.N; i++) {
    E = i * React_Xe_Xe.Estep + React_Xe_Xe.Emin;

    tmp = 0.;

    if (E / Ecoeff <= edata[0])
      tmp = Sdata[0];
    else if (E / Ecoeff >= edata[18])
      tmp = Sdata[18];
    else {
      int ind = 0;
      while (E / Ecoeff > edata[ind]) ind++;
      assert( ind >= 1 && "calculated index cases out of bounds access to Sdata" );
      tmp = Sdata[ind - 1] +
	    (Sdata[ind] - Sdata[ind - 1]) * (E / Ecoeff - edata[ind - 1]) /
		(edata[ind] - edata[ind - 1]);
    }

    React_Xe_Xe.CS[i] = CScoeff * tmp * 1.e4 * coll_fac_ntrl_ntrl * coll_fac_ntrl_ntrl;   // [1]

    fprintf(outputfile, "%f  %15.7e \n", E / Ecoeff, React_Xe_Xe.CS[i] / CScoeff);  // [eV] [cm^2] (scaled)
  }
  fclose(outputfile);


#if USE_TWOPLUS_IONS

  /******************************************************************************/
  /*                        Xe + e -> Xe++ + e + e + e                          */
  /*                         Ionization cross-section                           */
  /*             Source:  PHYSICAL REVIEW A VOLUME 35, NUMBER 3 p. 1041         */
  /******************************************************************************/
  // E <eV>           S <10^-16 cm2>  
  edata[0]  = 35;   Sdata[0]  = 0.013;
  edata[1]  = 36;   Sdata[1]  = 0.028;
  edata[2]  = 37;   Sdata[2]  = 0.045;
  edata[3]  = 38;   Sdata[3]  = 0.062;
  edata[4]  = 39;   Sdata[4]  = 0.079;
  edata[5]  = 40;   Sdata[5]  = 0.099;
  edata[6]  = 41;   Sdata[6]  = 0.116;
  edata[7]  = 42;   Sdata[7]  = 0.132;
  edata[8]  = 43;   Sdata[8]  = 0.149;
  edata[9]  = 44;   Sdata[9]  = 0.168;
  edata[10] = 45;  Sdata[10]  = 0.192;
  edata[11] = 46;  Sdata[11]  = 0.219;
  edata[12] = 47;  Sdata[12]  = 0.24 ;
  edata[13] = 48;  Sdata[13]  = 0.267;
  edata[14] = 49;  Sdata[14]  = 0.297;
  edata[15] = 50;  Sdata[15]  = 0.318;
  edata[16] = 51;  Sdata[16]  = 0.336;
  edata[17] = 52;  Sdata[17]  = 0.348;
  edata[18] = 53;  Sdata[18]  = 0.351;
  edata[19] = 54;  Sdata[19]  = 0.37 ;
  edata[20] = 55;  Sdata[20]  = 0.365;
  edata[21] = 56;  Sdata[21]  = 0.372;
  edata[22] = 57;  Sdata[22]  = 0.378;
  edata[23] = 58;  Sdata[23]  = 0.385;
  edata[24] = 59;  Sdata[24]  = 0.389;
  edata[25] = 60;  Sdata[25]  = 0.393;
  edata[26] = 61;  Sdata[26]  = 0.397;
  edata[27] = 62;  Sdata[27]  = 0.403;
  edata[28] = 63;  Sdata[28]  = 0.424;
  edata[29] = 64;  Sdata[29]  = 0.41 ;
  edata[30] = 65;  Sdata[30]  = 0.414;
  edata[31] = 66;  Sdata[31]  = 0.422;
  edata[32] = 67;  Sdata[32]  = 0.426;
  edata[33] = 68;  Sdata[33]  = 0.431;
  edata[34] = 69;  Sdata[34]  = 0.434;
  edata[35] = 70;  Sdata[35]  = 0.438;
  edata[36] = 71;  Sdata[36]  = 0.443;
  edata[37] = 72;  Sdata[37]  = 0.446;
  edata[38] = 73;  Sdata[38]  = 0.451;
  edata[39] = 74;  Sdata[39]  = 0.455;
  edata[40] = 75;  Sdata[40]  = 0.46;
  edata[41] = 76;  Sdata[41]  = 0.471;
  edata[42] = 77;  Sdata[42]  = 0.476;
  edata[43] = 78;  Sdata[43]  = 0.483;
  edata[44] = 79;  Sdata[44]  = 0.492;
  edata[45] = 80;  Sdata[45]  = 0.501;
  edata[46] = 81;  Sdata[46]  =  0.512;
  edata[47] = 82;  Sdata[47]  =  0.522;
  edata[48] = 83;  Sdata[48]  =  0.532;
  edata[49] = 84;  Sdata[49]  =  0.538;
  edata[50] = 85;  Sdata[50]  =  0.551;
  edata[51] = 86;  Sdata[51]  =  0.565;
  edata[52] = 87;  Sdata[52]  =  0.577;
  edata[53] =  88;  Sdata[53]  = 0.587; 
  edata[54] =  89;  Sdata[54]  = 0.592; 
  edata[55] =  90;  Sdata[55]  = 0.6; 
  edata[56] =  91;  Sdata[56]  = 0.609; 
  edata[57] =  92;  Sdata[57]  = 0.616; 
  edata[58] =  93;  Sdata[58]  = 0.621; 
  edata[59] =  94;  Sdata[59]  = 0.629; 
  edata[60] =  95;  Sdata[60]  = 0.634; 
  edata[61] =  96;  Sdata[61]  = 0.641; 
  edata[62] =  97;  Sdata[62]  = 0.647; 
  edata[63] =  98;  Sdata[63]  = 0.651; 
  edata[64] =  99;  Sdata[64]  = 0.646; 
  edata[65] = 100;  Sdata[65]  = 0.66; 

  Emin = 35;      //eV
  Eth =  35;      //eV
  Emax = 100.;    //eV

  Ecoeff = (dt/dr)*(dt/dr)*(1.+mn_over_me)/mn_over_me/(0.5*T_e0);
 
  CScoeff = collision_amplification_factor * ScaleF * n_e0 * dr_0 * ncoll_n / Ncell1;

  React_Xe_ipp.N = Nmax;
  React_Xe_ipp.Emin = Emin*Ecoeff;
  React_Xe_ipp.Eth =  Eth*Ecoeff;
  React_Xe_ipp.Emax = Emax*Ecoeff;
  React_Xe_ipp.Estep =  (React_Xe_ipp.Emax -  React_Xe_ipp.Emin)/React_Xe_ipp.N;

  outputfile = fopen("Xe_ipp.dat", "w");

  for ( int i = 0; i <= React_Xe_ipp.N; i++ )
  {
    E = i*React_Xe_ipp.Estep + React_Xe_ipp.Emin;

    tmp = 0.;

    if      (E/Ecoeff <= edata[0])  tmp = Sdata[0];  // TODO decide whether this behaviour is wanted
    else if (E/Ecoeff >= edata[65]) tmp = Sdata[65];
    else {
      int ind = 0;
      while (E/Ecoeff > edata[ind] ) ind++;
      tmp = Sdata[ind - 1] + (Sdata[ind]- Sdata[ind-1])*(E/Ecoeff - edata[ind-1])/(edata[ind]- edata[ind-1]);
    }

    React_Xe_ipp.CS[i] = CScoeff*tmp*1.e-16*coll_fac_ntrl_ntrl; // [1]

    fprintf(outputfile, "%f  %15.7e \n", E/Ecoeff, React_Xe_ipp.CS[i]/CScoeff ); // [eV]  [cm^2] (scaled)
  }
  fclose(outputfile); 


  /******************************************************************************/
  /*                        Xe+ + e -> Xe++ + e + e                             */
  /*                         Ionization cross-section                           */
  /*             Source:  Presentation of Alfred Müller (2009)    ???           */
  /******************************************************************************/
  // ???       E <eV>                S <10^-22 m2>      ????
	edata[0]=21.042 ;	  Sdata[0]=11.8354;
	edata[1]=22.0825;	  Sdata[1]=29.732;
	edata[2]=23.0799;	  Sdata[2]=63.7726;
	edata[3]=24.1134;	  Sdata[3]= 88.6765;
	edata[4]=25.0922;	  Sdata[4]= 115.752;
	edata[5]=26.1093;	  Sdata[5]= 141.201;
	edata[6]=27.1126;	  Sdata[6]= 162.425;
	edata[7]=28.155 ;	  Sdata[7]= 191.97;
	edata[8]=29.1195;	  Sdata[8]= 204.51;
	edata[9]=30.1768;	  Sdata[9]= 218.855;
	edata[10]=31.1493;	  Sdata[10]= 236.329;
	edata[11]=32.1523;	  Sdata[11]= 248.377;
	edata[12]=33.1219;	  Sdata[12]= 259.864;
	edata[13]=34.1196;	  Sdata[13]= 261.055;
	edata[14]=35.0784;	  Sdata[14]= 267.029;
	edata[15]=36.2066;	  Sdata[15]= 268.254;
	edata[16]=37.0759;	  Sdata[16]= 264.658;
	edata[17]=38.1167;	  Sdata[17]= 262.293;
	edata[18]=39.1101;	  Sdata[18]= 265.883;
	edata[19]=40.1285;	  Sdata[19]= 262.319;
	edata[20]=42.2457;	  Sdata[20]= 257.65 ;
	edata[21]=44.2107;	  Sdata[21]= 247.413;
	edata[22]=46.1779;	  Sdata[22]= 251.943;
	edata[23]=48.2314;	  Sdata[23]= 248.575;
	edata[24]=50.2779;	  Sdata[24]= 253.126;
	edata[25]=55.2859;	  Sdata[25]= 253.173;
	edata[26]=60.3122;	  Sdata[26]= 246.449;
	edata[27]=65.4089;	  Sdata[27]= 252.115;
	edata[28]=70.374 ;	  Sdata[28]= 244.307;
	edata[29]=75.2712;	  Sdata[29]= 248.792;
	edata[30]=80.5066;	  Sdata[30]= 244.372;
	edata[31]=85.4302;	  Sdata[31]= 248.854;
	edata[32]=90.6518;	  Sdata[32]= 243.327;
	edata[33]=95.4359;	  Sdata[33]= 242.256;
	edata[34]=100.471;	  Sdata[34]= 236.872;
	edata[35]=110.694;	  Sdata[35]= 229.547;
	edata[36]=120.758;	  Sdata[36]= 223.45 ;
	edata[37]=130.696;	  Sdata[37]= 212.657;
	edata[38]=140.615;	  Sdata[38]= 204.218;
	edata[39]=150.99 ;	  Sdata[39]= 197.892;
	edata[40]=160.855;	  Sdata[40]= 194.375;
	edata[41]=171.02 ;	  Sdata[41]= 184.146;
	edata[42]=181.47 ;	  Sdata[42]= 177.634;
	edata[43]=191.045;	  Sdata[43]= 174.473;
	edata[44]=200.724;	  Sdata[44]= 168.3  ;
	edata[45]=221.148;	  Sdata[45]= 162.36 ;
	edata[46]=241.724;	  Sdata[46]= 153.132;
	edata[47]=261.098;	  Sdata[47]= 144.424;
	edata[48]=281.47 ;	  Sdata[48]= 137.447;
	edata[49]=301.638;	  Sdata[49]= 132.589;
	edata[50]=331.667;	  Sdata[50]= 125.055;
	edata[51]=361.809;	  Sdata[51]= 117.415;
	edata[52]=402.571;	  Sdata[52]= 108.763;
	edata[53]=451.477;	  Sdata[53]= 98.9467;
	edata[54]=502.347;	  Sdata[54]= 92.9054;
	edata[55]=552.343;	  Sdata[55]= 85.2839;
	edata[56]=602.54 ;	  Sdata[56]= 80.0738;
	edata[57]=653.412;	  Sdata[57]= 75.5212;
	edata[58]=702.99 ;	  Sdata[58]= 70.5862;
	edata[59]=751.883;	  Sdata[59]= 69.0197;
	edata[60]=804.158;	  Sdata[60]= 65.6843;
	edata[61]=901.857;	  Sdata[61]= 60.2982;
	edata[62]=1005.45;	  Sdata[62]= 55.6036;
					   
  Emin = edata[0];
  Eth = edata[0];
  Emax = edata[62];
 
  Ecoeff = (dt/dr)*(dt/dr)*(1.+mi_over_me)/mi_over_me/(0.5*T_e0);

  CScoeff = collision_amplification_factor * ScaleF * n_e0 * dr_0 * ncoll_n / Ncell1; 
  
  React_Xep_ipp.N = Nmax;
  React_Xep_ipp.Emin = Emin*Ecoeff;
  React_Xep_ipp.Eth =  Eth*Ecoeff;
  React_Xep_ipp.Emax = Emax*Ecoeff;
  React_Xep_ipp.Estep =  (React_Xep_ipp.Emax -  React_Xep_ipp.Emin)/React_Xep_ipp.N;

  outputfile = fopen("Xep_ipp.dat", "w");

  for ( int i = 0; i <= React_Xep_ipp.N; i++ )
  {
    E = i*React_Xep_ipp.Estep + React_Xep_ipp.Emin;

    tmp = 0.;

    if      (E/Ecoeff <= Emin)  tmp = Sdata[0]; // TODO decide whether this behaviour is wanted
    else if (E/Ecoeff >= Emax) tmp = Sdata[62];
    else  {
      int ind = 0;
      while (E/Ecoeff > edata[ind] ) ind++;
      tmp = Sdata[ind - 1] + (Sdata[ind]- Sdata[ind-1])*(E/Ecoeff - edata[ind-1])/(edata[ind]- edata[ind-1]);
    }

    React_Xep_ipp.CS[i] = CScoeff*tmp*1.e-18;

    fprintf(outputfile, "%f  %15.7e \n", E/Ecoeff, React_Xep_ipp.CS[i]/CScoeff );  // should be [eV] [cm^2] -> data in [eV] [10^-22 m^2] ??
  }
  fclose(outputfile); 

#endif

}//end: init_reactions_Xe

/*******************************************************************************
		   Here we init records with reactions information for Argon
*******************************************************************************/
void init_reactions_Ar(void) {
  double E;
  int i, j;
  double tmp;
  double Ecoeff, CScoeff;
  double Emin, Emax, Eth;

  FILE *outputfile;

  double edata[80], Sdata[80];

  int Nmax;

  Nmax = Npoints_CS - 1;


  /********************************************************************************/
  /*                            Ar + e -> Ar+ + 2e */
  /*                         Ionization cross-section */
  /*                    Source:  Konstantin 1D PIC code */
  /********************************************************************************/

  // use function: coll_el_ntrl_ioniz()
  // neutrals are scaled onto electron step, no additional scaling needed

  Emin = 15.8;  // eV
  Eth = 15.8;  // eV
  Emax = 50.;  // eV

  Ecoeff = (dt / dr) * (dt / dr) * (1. + mn_over_me) / mn_over_me / (0.5 * T_e0);

  CScoeff = collision_amplification_factor * ScaleF * n_e0 * dr_0 * ncoll_n / Ncell1; 
  
  React_Ar_i.N = Nmax;
  React_Ar_i.Emin = Emin * Ecoeff;
  React_Ar_i.Eth = Eth * Ecoeff;
  React_Ar_i.Emax = Emax * Ecoeff;
  React_Ar_i.Estep = (React_Ar_i.Emax - React_Ar_i.Emin) / React_Ar_i.N;

  outputfile = fopen("Ar_i.dat", "w");

  for (i = 0; i <= React_Ar_i.N; i++) {
    E = i * React_Ar_i.Estep + React_Ar_i.Emin;

	tmp = 9.7E-14 * (E / Ecoeff - 15.8) / pow((70 + E / Ecoeff), 2) +
	6.E-22 * pow((E / Ecoeff - 15.8), 2) * exp(-E / Ecoeff / 9);
    
	React_Ar_i.CS[i] = CScoeff * tmp * coll_fac_ntrl_ntrl;


    fprintf(outputfile, "%f  %15.7e \n", E / Ecoeff,
	    React_Ar_i.CS[i] / CScoeff);
  }
  fclose(outputfile);

  /********************************************************************************/
  /*                             Ar + e -> Ar + e */
  /*                         Total elastic cross-section */
  /*                     Source:  Konstantin 1D PIC code */
  /********************************************************************************/
  
  //    E <eV>                  S <10^-16 cm2>
  edata[0] = 0.0010;	      Sdata[0] = 7.50;
  edata[1] = 0.0020;	      Sdata[1] = 7.10;
  edata[2] = 0.0030;	      Sdata[2] = 6.70;
  edata[3] = 0.0050;	      Sdata[3] = 6.10;
  edata[4] = 0.0070;	      Sdata[4] = 5.40;
  edata[5] = 0.0085;	      Sdata[5] = 5.05;
  edata[6] = 0.0100;	      Sdata[6] = 4.60;
  edata[7] = 0.0150;	      Sdata[7] = 3.75;
  edata[8] = 0.0200;	      Sdata[8] = 3.25;
  edata[9] = 0.0300;	      Sdata[9] = 2.50;
  edata[10] = 0.0400;	      Sdata[10] = 2.05;
  edata[11] = 0.0500;	      Sdata[11] = 1.73;
  edata[12] = 0.0700;	      Sdata[12] = 1.130;
  edata[13] = 0.1000;         Sdata[13] = 0.590;
  edata[14] = 0.1200;	      Sdata[14] = 0.400;
  edata[15] = 0.1500;	      Sdata[15] = 0.230;
  edata[16] = 0.1700;	      Sdata[16] = 0.160;
  edata[17] = 0.2000;	      Sdata[17] = 0.103;
  edata[18] = 0.2500;	      Sdata[18] = 0.091;
  edata[19] = 0.3000;	      Sdata[19] = 0.153;
  edata[20] = 0.3500;	      Sdata[20] = 0.235;
  edata[21] = 0.4000;	      Sdata[21] = 0.33;
  edata[22] = 0.5000;	      Sdata[22] = 0.51;
  edata[23] = 0.7000;	      Sdata[23] = 0.86;
  edata[24] = 1.0000;	      Sdata[24] = 1.38;
  edata[25] = 1.2000;	      Sdata[25] = 1.66;
  edata[26] = 1.3000;	      Sdata[26] = 1.82;
  edata[27] = 1.5000;	      Sdata[27] = 2.10;
  edata[28] = 1.7000;	      Sdata[28] = 2.3;
  edata[29] = 1.9000;	      Sdata[29] = 2.5;
  edata[30] = 2.1000;	      Sdata[30] = 2.8;
  edata[31] = 2.2000;	      Sdata[31] = 2.9;
  edata[32] = 2.5000;	      Sdata[32] = 3.3;
  edata[33] = 2.8000;	      Sdata[33] = 3.8;
  edata[34] = 3.0000;	      Sdata[34] = 4.1;
  edata[35] = 3.3000;	      Sdata[35] = 4.5;
  edata[36] = 3.6000;	      Sdata[36] = 4.9;
  edata[37] = 4.0000;	      Sdata[37] = 5.4;
  edata[38] = 4.5000;	      Sdata[38] = 6.1;
  edata[39] = 5.0000;	      Sdata[39] = 6.7;
  edata[40] = 6.0000;	      Sdata[40] = 8.1;
  edata[41] = 7.0000;	      Sdata[41] = 9.6;
  edata[42] = 8.0000;	      Sdata[42] = 11.7;
  edata[43] = 10.0000;	      Sdata[43] = 15.0;
  edata[44] = 12.0000;	      Sdata[44] = 15.2;
  edata[45] = 15.0000;	      Sdata[45] = 14.1;
  edata[46] = 17.0000;	      Sdata[46] = 13.1;
  edata[47] = 20.0000;	      Sdata[47] = 11.0;
  edata[48] = 25.0000;	      Sdata[48] = 9.45;
  edata[49] = 30.0000;	      Sdata[49] = 8.74;
  edata[50] = 50.0000;	      Sdata[50] = 6.90;
  edata[51] = 75.0000;	      Sdata[51] = 5.85;
  edata[52] = 100.0000;	      Sdata[52] = 5.25;
  edata[53] = 150.0000;	      Sdata[53] = 4.24;
  edata[54] = 200.0000;	      Sdata[54] = 3.76;
  edata[55] = 300.0000;	      Sdata[55] = 3.02;

  // use function: coll_el_all_fake() 
  // second input partner (should be electrons, m_e is hardwired) loses some energy
  // neutrals are scaled onto electron step, no additional scaling needed
  
  Emin = 0.0;  // eV
  Eth = 0.0;  // eV
  Emax = 100.;  // eV

  Ecoeff =
      (dt / dr) * (dt / dr) * (1. + mn_over_me) / mn_over_me / (0.5 * T_e0);

  CScoeff = collision_amplification_factor * ScaleF * n_e0 * dr_0 * ncoll_n / Ncell1; 
  
  React_Ar_el.N = Nmax;
  React_Ar_el.Emin = Emin * Ecoeff;
  React_Ar_el.Eth = Eth * Ecoeff;
  React_Ar_el.Emax = Emax * Ecoeff;
  React_Ar_el.Estep = (React_Ar_el.Emax - React_Ar_el.Emin) / React_Ar_el.N;

  outputfile = fopen("Ar_el.dat", "w");

  for (i = 0; i <= React_Ar_el.N; i++) {
    E = i * React_Ar_el.Estep + React_Ar_el.Emin;
    tmp = 0.;

    if (E / Ecoeff <= edata[0])
      tmp = Sdata[0];
    else if (E / Ecoeff >= edata[55])
      tmp = Sdata[55];
    else {
      int ind = 0;
      while (E / Ecoeff > edata[ind]) ind++;
      assert( ind >= 1 && "calculated index cases out of bounds access to Sdata" );
      tmp = Sdata[ind - 1] +
	    (Sdata[ind] - Sdata[ind - 1]) * (E / Ecoeff - edata[ind - 1]) /
		(edata[ind] - edata[ind - 1]);
    }

    React_Ar_el.CS[i] = CScoeff * tmp * 1.e-16 * coll_fac_ntrl_ntrl;

    fprintf(outputfile, "%f  %15.7e \n", E / Ecoeff,
	    React_Ar_el.CS[i] / CScoeff);
  }
  fclose(outputfile);

  /********************************************************************************/
  /*                      Ar + e -> Ar* + e */
  /*                Total excitation cross-section */
  /*       Source:  Konstantin 1D PIC code */
  /********************************************************************************/

  edata[0] = 0.0000;		Sdata[0] = 0.0000;
  edata[1] = 11.5000;		Sdata[1] = 0.0000;
  edata[2] = 12.7000;		Sdata[2] = 0.0700;
  edata[3] = 13.7000;		Sdata[3] = 0.1410;
  edata[4] = 14.7000;		Sdata[4] = 0.2280;
  edata[5] = 15.9000;		Sdata[5] = 0.3800;
  edata[6] = 16.5000;		Sdata[6] = 0.4800;
  edata[7] = 17.5000;		Sdata[7] = 0.6100;
  edata[8] = 18.5000;		Sdata[8] = 0.7500;
  edata[9] = 19.9000;		Sdata[9] = 0.9200;
  edata[10] = 22.2000;		Sdata[10] = 1.1700;
  edata[11] = 24.7000;		Sdata[11] = 1.3300;
  edata[12] = 27.0000;		Sdata[12] = 1.4200;
  edata[13] = 30.0000;		Sdata[13] = 1.4400;
  edata[14] = 33.0000;		Sdata[14] = 1.4100;
  edata[15] = 35.3000;		Sdata[15] = 1.3400;
  edata[16] = 42.0000;		Sdata[16] = 1.2500;
  edata[17] = 48.0000;		Sdata[17] = 1.1600;
  edata[18] = 52.0000;		Sdata[18] = 1.1100;
  edata[19] = 70.;		Sdata[19] = 0.94;
  edata[20] = 100.;		Sdata[20] = 0.76;
  edata[21] = 150.;		Sdata[21] = 0.60;
  edata[22] = 200.;		Sdata[22] = 0.505;
  edata[23] = 300.;		Sdata[23] = 0.395;
  edata[24] = 500.;		Sdata[24] = 0.28;
  edata[25] = 700.;		Sdata[25] = 0.225;
  edata[26] = 1000.;		Sdata[26] = 0.177;
  edata[27] = 1500.;		Sdata[27] = 0.136;
  edata[28] = 2000.;		Sdata[28] = 0.11;
  edata[29] = 3000.;		Sdata[29] = 0.083;
  edata[30] = 5000.;		Sdata[30] = 0.058;
  edata[31] = 7000.;		Sdata[31] = 0.045;
  edata[32] = 10000.;		Sdata[32] = 0.035;

  // use function: coll_el_all_fake() 
  // second input partner (should be electrons, m_e is hardwired) loses some energy
  // neutrals are scaled onto electron step, no additional scaling needed
  
  Emin = 11.5;  // eV
  Eth = 11.5;  // eV
  Emax = 7000.;  // eV

  Ecoeff =
      (dt / dr) * (dt / dr) * (1. + mn_over_me) / mn_over_me / (0.5 * T_e0);

  CScoeff = collision_amplification_factor * ScaleF * n_e0 * dr_0 * ncoll_n / Ncell1; 
  
  React_Ar_tex.N = Nmax;
  React_Ar_tex.Emin = Emin * Ecoeff;
  React_Ar_tex.Eth = Eth * Ecoeff;
  React_Ar_tex.Emax = Emax * Ecoeff;
  React_Ar_tex.Estep = (React_Ar_tex.Emax - React_Ar_tex.Emin) / React_Ar_tex.N;

  outputfile = fopen("Ar_tex.dat", "w");

  for (i = 0; i <= React_Ar_tex.N; i++) {
    E = i * React_Ar_tex.Estep + React_Ar_tex.Emin;

    tmp = 0.;

    if (E / Ecoeff <= edata[1])
      tmp = Sdata[1];
    else if (E / Ecoeff >= edata[32])
      tmp = Sdata[32];
    else {
      int ind = 0;
      while (E / Ecoeff > edata[ind]) ind++;
      assert( ind >= 1 && "calculated index cases out of bounds access to Sdata" );
      tmp = Sdata[ind - 1] +
	    (Sdata[ind] - Sdata[ind - 1]) * (E / Ecoeff - edata[ind - 1]) /
		(edata[ind] - edata[ind - 1]);
    }

    React_Ar_tex.CS[i] = CScoeff * tmp * 1.e-16 * coll_fac_ntrl_ntrl;

    fprintf(outputfile, "%f  %15.7e \n", E / Ecoeff,
	    React_Ar_tex.CS[i] / CScoeff);
  }
  fclose(outputfile);

  /********************************************************************************/
  /*                      Ar + Ar+ -> Ar + Ar+ */
  /*        Elastic and Charge Exchange  (isotropic scattering aproximatiuon) */
  /* Source: http://jila.colorado.edu/~avp/collision_data/ionneutral/IONATOM.TXT
  */
  /********************************************************************************/

  // use function: coll_ion_ntrl()
  // Ion is scaled to dt_ntrl, neutral remains in its velocity scale
  // Reactionrate has to compensate higher velocity: CScoeff0 / dt_ntrl
  // Energy is compared at neutral scale: Ecoeff0 * dt_ntrl * dt_ntrl

  Emin = 0.0;  // eV
  Eth = 0.0;  // eV
  Emax = 100.0;  // eV

  Ecoeff = (dt / dr) * (dt / dr) * (mi_over_me + mn_over_me) / mn_over_me /
	   mi_over_me / (0.5 * T_e0) * dt_ntrl *
	   dt_ntrl;  // we bring ion velocity on the neutral scale

  CScoeff = collision_amplification_factor * ScaleF * n_e0 * dr_0 * ncoll_n / Ncell1 /
	    dt_ntrl;  // we bring ion velocity on the neutral scale

  React_Arp_Ar.N = Nmax;
  React_Arp_Ar.Emin = Emin * Ecoeff;
  React_Arp_Ar.Eth = Eth * Ecoeff;
  React_Arp_Ar.Emax = Emax * Ecoeff;
  React_Arp_Ar.Estep = (React_Arp_Ar.Emax - React_Arp_Ar.Emin) / React_Arp_Ar.N;

  outputfile = fopen("Arp_Ar.dat", "w");

  for (i = 0; i <= React_Arp_Ar.N; i++) {
    E = i * React_Arp_Ar.Estep + React_Arp_Ar.Emin;

    if (E == 0.) E = 0.000092;  // to avoid inf at E=0

    tmp = 1.15E-14 * pow((E / Ecoeff), -0.1) * pow((1 + 0.015 * Ecoeff / E), 0.6);

    React_Arp_Ar.CS[i] = CScoeff * tmp * coll_fac_ntrl_ntrl;

    fprintf(outputfile, "%f  %15.7e \n", E / Ecoeff,
	    React_Arp_Ar.CS[i] / CScoeff);
  }
  fclose(outputfile);

  /***************************************************************************************************/
  /*                      Ar + Ar ->  Ar + Ar */
  /*                Atom-Atom elastic collision */
  /*   Source:
   * http://jila.colorado.edu/~avp/collision_data/neutralneutral/atomatom.txt */
  /*******************3********************************************************************************/

  edata[0] = 1.00000E-02;	  Sdata[0] = 4.15120E-18;
  edata[1] = 1.77828E-02;	  Sdata[1] = 3.89469E-18;
  edata[2] = 3.16228E-02;	  Sdata[2] = 3.58990E-18;
  edata[3] = 5.62341E-02;	  Sdata[3] = 2.88137E-18;
  edata[4] = 1.00000E-01;	  Sdata[4] = 2.95445E-18;
  edata[5] = 1.77828E-01;	  Sdata[5] = 2.39723E-18;
  edata[6] = 3.16228E-01;	  Sdata[6] = 2.10701E-18;
  edata[7] = 5.62341E-01;	  Sdata[7] = 2.21634E-18;
  edata[8] = 1.00000E+00;	  Sdata[8] = 2.06321E-18;
  edata[9] = 1.77828E+00;	  Sdata[9] = 1.69308E-18;
  edata[10] = 3.16228E+00;	  Sdata[10] = 1.32380E-18;
  edata[11] = 5.62341E+00;	  Sdata[11] = 1.04412E-18;
  edata[12] = 1.00000E+01;	  Sdata[12] = 8.58338E-19;
  edata[13] = 1.77828E+01;	  Sdata[13] = 7.41028E-19;
  edata[14] = 3.16228E+01;	  Sdata[14] = 6.67289E-19;
  edata[15] = 5.62341E+01;	  Sdata[15] = 6.18802E-19;
  edata[16] = 1.00000E+02;	  Sdata[16] = 5.84667E-19;
  edata[17] = 1.77828E+02;	  Sdata[17] = 5.58245E-19;
  edata[18] = 3.16228E+02;	  Sdata[18] = 5.35810E-19;
  edata[19] = 5.62341E+02;	  Sdata[19] = 5.14848E-19;
  edata[20] = 1.00000E+03;	  Sdata[20] = 4.94013E-19;
  edata[21] = 1.77828E+03;	  Sdata[21] = 4.72663E-19;
  edata[22] = 3.16228E+03;	  Sdata[22] = 4.50369E-19;
  edata[23] = 5.62341E+03;	  Sdata[23] = 4.26661E-19;
  edata[24] = 1.00000E+04;	  Sdata[23] = 4.02005E-19;

 // use function: coll_ntrl_ntrl()
 // neutrals are kept in their scale
 // Reactionrate has to compensate higher velocity: CScoeff0 / dt_ntrl
 // Energy is compared at neutral scale: Ecoeff0 * dt_ntrl * dt_ntrl
 // Threshold is used as a maximum Sigma_v, thus has to be scaled 
  
  Emin = 1.0E-2;  // eV
  Emax = 3.16228E+02;  // eV

  CScoeff = collision_amplification_factor * ScaleF * n_e0 * dr_0 * ncoll_n_n / Ncell1 / dt_ntrl; 

  Ecoeff = (dt / dr) * (dt / dr) * 2. / mn_over_me / (0.5 * T_e0) * dt_ntrl *
	   dt_ntrl;

  React_Ar_Ar.N = Nmax;
  React_Ar_Ar.Emin = Emin * Ecoeff;

  React_Ar_Ar.Emax = Emax * Ecoeff;
  React_Ar_Ar.Estep = (React_Ar_Ar.Emax - React_Ar_Ar.Emin) / React_Ar_Ar.N;

  React_Ar_Ar.Eth = sqrt(edata[18] * Ecoeff) * Sdata[18] * 1.e4 * CScoeff * coll_fac_ntrl_ntrl * coll_fac_ntrl_ntrl;

 //  printf("-> Ar-Ar MAX(Sigma*v) %f \n",  React_Ar_Ar.Eth );

  outputfile = fopen("Ar_Ar.dat", "w");

  for (i = 0; i <= React_Ar_Ar.N; i++) {
    E = i * React_Ar_Ar.Estep + React_Ar_Ar.Emin;

    tmp = 0.;

    if (E / Ecoeff <= edata[0])
      tmp = Sdata[0];
    else if (E / Ecoeff >= edata[18])
      tmp = Sdata[18];
    else {
      int ind = 0;
      while (E / Ecoeff > edata[ind]) ind++;
      assert( ind >= 1 && "calculated index cases out of bounds access to Sdata" );
      tmp = Sdata[ind - 1] +
	    (Sdata[ind] - Sdata[ind - 1]) * (E / Ecoeff - edata[ind - 1]) /
		(edata[ind] - edata[ind - 1]);
    }

    React_Ar_Ar.CS[i] = CScoeff * tmp * 1.e4 * coll_fac_ntrl_ntrl * coll_fac_ntrl_ntrl;  // karl 22.10. additional collision factor

    fprintf(outputfile, "%f  %15.7e \n", E / Ecoeff,
	    React_Ar_Ar.CS[i] / CScoeff);
  }
  fclose(outputfile);
}


/*******************************************************************************
		   Here we init records with reactions information for Oxygen
*******************************************************************************/
void init_reactions_O2(void)
{
    double E, Ldb, Omega_pe;
    int i, j;
    double tmp;
    double Ecoeff, CScoeff;
    double Emin, Emax, Eth, Efc;
    double n_pressureX, ampl_pressure;
    FILE *outputfile;

  double edata[80], Sdata[80];

  int Npoints;

  Npoints=Npoints_CS-1;

   Ldb = 7.43e2*sqrt(T_e0/n_e0);    //T_e = 1
   printf("Debye-length= %f \n", Ldb);

   Omega_pe = 5.64e4*sqrt(n_e0);
 
   CScoeff = collision_amplification_factor * ScaleF* n_e0*dr_0*ncoll_n/Ncell1;      
#if USE_FEEDGAS_OXYGEN // this should be an own RF flag 
   n_pressureX= pressure  * 6.2418e12/Ti_over_Te; // neutrals per cm³ (ideal gas law)
   ampl_pressure = n_pressureX/n_e0; // amplification factor
   ampl_pressure = 1;                //TODO check why it is set 1 here, nullifies pressure set
#else
   ampl_pressure = 1;
#endif
 printf("ampl_pressure: %f\n",ampl_pressure); 
  /*************************************************************************************************
    ellastic scattering of e- on O2
  ***************************************************************************************************/
   
  edata[0]   =   0.;                  Sdata[0]   = 0.35; 
  edata[1]   =   1.e-3;             	Sdata[1]   = 0.35;
  edata[2]   =   0.002;            	Sdata[2]   = 0.36;  
  edata[3]   =   0.003;            	Sdata[3]   = 0.4;     
  edata[4]   =   0.005;            	Sdata[4]   = 0.5;     
  edata[5]   =   0.007;            	Sdata[5]   = 0.58;   
  edata[6]   =   0.0085;          	Sdata[6]   = 0.64;     
  edata[7]   =   0.01;              	Sdata[7]   = 0.7;     
  edata[8]   =   0.015;             	Sdata[8]   = 0.87;   
  edata[9]   =   0.02;               	Sdata[9]   = 0.99;  
  edata[10]  =  0.03;               	Sdata[10]  =1.24;    
  edata[11]  =  0.04;               	Sdata[11]  =1.44;  
  edata[12]  =  0.05;               	Sdata[12]  =1.6;        
  edata[13]  =  0.07;              	Sdata[13]  =2.1;       
  edata[14]  =  0.1;                	Sdata[14]  =2.5;      
  edata[15]  =  0.12;              	Sdata[15]  =2.8;      
  edata[16]  =  0.15;              	Sdata[16]  =3.1;      
  edata[17]  =  0.17;              	Sdata[17]  =3.3;      
  edata[18]  =  0.2;                	Sdata[18]  =3.6;    
  edata[19]  =  0.25;              	Sdata[19]  =4.1;        
  edata[20]  =  0.3;                	Sdata[20]  =4.5;        
  edata[21]  =  0.35;              	Sdata[21]  =4.7;      
  edata[22]  =  0.4;                	Sdata[22]  =5.2;      
  edata[23]  =  0.5;               	Sdata[23]  =5.7;       
  edata[24]  =  0.7;                	Sdata[24]  =6.1;    
  edata[25]  =  1.;                   Sdata[25]  =7.2;      
  edata[26]  =  1.2;                 	Sdata[26]  =7.9;    
  edata[27]  =  1.3;                	Sdata[27]  =7.9;      
  edata[28]  =  1.5;                	Sdata[28]  =7.6;      
  edata[29]  =  1.7;                 	Sdata[29]  =7.3;    
  edata[30]  =  1.9;                	Sdata[30]  =6.9;      
  edata[31]  =  2.1;                 	Sdata[31]  =6.6;      
  edata[32]  =  2.2;       			Sdata[32]  =6.5;      
  edata[33]  =  2.5;                 	Sdata[33]  =6.1;    
  edata[34]  =  2.8;                 	Sdata[34]  =5.8;      
  edata[35]  =  3.;                 	Sdata[35]  =5.7;    
  edata[36]  =  3.3;                	Sdata[36]  =5.5;     
  edata[37]  =  3.6;                	Sdata[37]  =5.45;     
  edata[38]  =  4.;                   Sdata[38]  =5.5;      
  edata[39]  =  4.5;                	Sdata[39]  =5.55;    
  edata[40]  =  5.;                   Sdata[40]  =5.6;         
  edata[41]  =  6.;                   Sdata[41]  =6.;                              
  edata[42]  =  7.;                   Sdata[42]  =6.6;                          
  edata[43]  =  8.;                  	Sdata[43]  =7.1;          
  edata[44]  =  10.;                	Sdata[44]  =8.;             
  edata[45]  =  12.;                	Sdata[45]  =8.5;                                        
  edata[46]  =  15.;                	Sdata[46]  =8.8;      
  edata[47]  =  17.;               	Sdata[47]  =8.7;             
  edata[48]  =  20.;               	Sdata[48]  =8.6;             
  edata[49]  =  25.;               	Sdata[49]  =8.2;         
  edata[50]  =  30.;               	Sdata[50]  =8.; 
	     
  edata[51]  =  50.;               	Sdata[51]  =7.7;         
  edata[52]  =  75.;               	Sdata[52]  =6.8;             
  edata[53]  =  100.;             	Sdata[53]  =6.5;             
  edata[54]  =  150.;             	Sdata[54]  =6.7;      


   
   Emin = 0.;             //eV
   Eth =  0.;            //eV
   Emax = 30.;         //eV
   Efc =  0.;          //eV

   Ecoeff = (dt/dr)*(dt/dr)*(1.+mn_over_me)/mn_over_me/(0.5*T_e0);

   React_O2_el.Emin = Emin*Ecoeff;
   React_O2_el.Eth =  Eth*Ecoeff;
   React_O2_el.Emax = Emax*Ecoeff;
  // React_O2_el.Wfc =  sqrt(Efc*Ecoeff);
   React_O2_el.N = Npoints;
   React_O2_el.Estep =  (React_O2_el.Emax -  React_O2_el.Emin)/Npoints;

    outputfile = fopen("xsct/O2_e_el.dat", "w");

   for ( i = 0; i <= React_O2_el.N; i++ )
   {
   E = i*React_O2_el.Estep + React_O2_el.Emin;

   tmp = 0.;

   if      (E/Ecoeff <= edata[0])  tmp = Sdata[0];
   else if (E/Ecoeff >= edata[54]) tmp = Sdata[54];
       else
     {
       int ind = 0;
      while (E/Ecoeff > edata[ind] ) ind++;
      assert( ind >= 1 && "calculated index cases out of bounds access to Sdata" );
      tmp = Sdata[ind - 1] + (Sdata[ind]- Sdata[ind-1])*(E/Ecoeff - edata[ind-1])/(edata[ind]- edata[ind-1]);
     }

   React_O2_el.CS[i] = CScoeff*tmp*1.e-16*coll_fac_ntrl_ntrl*ampl_pressure;     

      fprintf(outputfile, "%f  %15.7e \n", E/Ecoeff, React_O2_el.CS[i]/CScoeff );
     }
   fclose(outputfile); 
   
  /********************************************************************************************* 
   e- + O2 -> O2+  + 2e-
 *********************************************************************************************/

  edata[0]   = 12.06; 	    Sdata[0]   = 0.;        
  edata[1]   = 13.;       	Sdata[1]   = 0.023;
  edata[2]   = 18.;       	Sdata[2]   = 0.2;      
  edata[3]   = 28.;       	Sdata[3]   = 0.74;   
  edata[4]   = 38.;       	Sdata[4]   = 1.32;   
  edata[5]   = 48.;       	Sdata[5]   = 1.8;     
  edata[6]   = 58.;       	Sdata[6]   = 2.1;     
  edata[7]   = 68.;       	Sdata[7]   = 2.33;   
  edata[8]   = 78.;        	Sdata[8]   = 2.5;     
  edata[9]   = 88.;       	Sdata[9]   = 2.6;      
  edata[10]  = 100.;      	Sdata[10]  = 2.7;      
  edata[11]  = 150.;      	Sdata[11]  = 2.7;      



   Emin = 12.06;             //eV
   Eth =  12.06;            //eV
   Emax = 150.;         //eV
   Efc =  0.;          //eV

   Ecoeff = (dt/dr)*(dt/dr)*(1.+mn_over_me)/mn_over_me/(0.5*T_e0);

   React_O2_i.Emin = Emin*Ecoeff;
 React_O2_i.Eth =  Eth*Ecoeff;
 React_O2_i.Emax = Emax*Ecoeff;
 //React_O2_i.Wfc =  sqrt(Efc*Ecoeff);
 React_O2_i.N = Npoints;
 React_O2_i.Estep =  (React_O2_i.Emax -  React_O2_i.Emin)/Npoints;

  outputfile = fopen("xsct/O2_i.dat", "w");

 for ( i = 0; i <= React_O2_i.N; i++ )
 {
 E = i*React_O2_i.Estep + React_O2_i.Emin;

 tmp = 0.;

 if      (E/Ecoeff <= edata[0])  tmp = Sdata[0];
 else if (E/Ecoeff >= edata[11]) tmp = Sdata[11];
     else
   {
     int ind = 0;
    while (E/Ecoeff > edata[ind] ) ind++;
    assert( ind >= 1 && "calculated index cases out of bounds access to Sdata" );
    tmp = Sdata[ind - 1] + (Sdata[ind]- Sdata[ind-1])*(E/Ecoeff - edata[ind-1])/(edata[ind]- edata[ind-1]);
   }

 React_O2_i.CS[i] = CScoeff*tmp*1.e-16*coll_fac_ntrl_ntrl*ampl_pressure;

  fprintf(outputfile, "%f  %15.7e \n", E/Ecoeff, React_O2_i.CS[i]/CScoeff/coll_fac_ntrl_ntrl/ampl_pressure);
   }
 fclose(outputfile); 
 
/****************************************************************************************
  e + O2 -> 2O + e (neutral dissociation)
****************************************************************************************/

edata[0]   = 8.;         	Sdata[0]   = 0.;       
edata[1]   = 13.5;     		Sdata[1]   = 0.22; 
edata[2]   = 18.5;     		Sdata[2]   = 0.53;   
edata[3]   = 21.;        	Sdata[3]   = 0.56;  
edata[4]   = 23.5;     		Sdata[4]   = 0.52;  
edata[5]   = 28.5;     		Sdata[5]   = 0.59;  
edata[6]   = 33.5;     		Sdata[6]   = 0.66;  
edata[7]   = 38.5;     		Sdata[7]   = 0.61;  
edata[8]   = 48.5;    		Sdata[8]   = 0.53;  
edata[9]   = 58.5;    		Sdata[9]   = 0.44;   
edata[10]  = 73.5;      	Sdata[10]  = 0.37;   
edata[11]  = 98.5;      	Sdata[11]  = 0.33;   
edata[12]  = 148.5;    		Sdata[12]  = 0.3;     



 Emin = 8.;             //eV
 Eth =  8.;            //eV
 Emax = 148.5;         //eV
 Efc =  1.5;          //eV


//Ecoeff = (Omega_pe/dx)*(Omega_pe/dx)*(1.+M_ions[13])/M_ions[13]/(0.5*T_e);

 Ecoeff = (dt/dr)*(dt/dr)*(1.+mn_over_me)/mn_over_me/(0.5*T_e0);

 React_O2_d.Emin = Emin*Ecoeff;
 React_O2_d.Eth =  Eth*Ecoeff;
 React_O2_d.Emax = Emax*Ecoeff;
// React_O2_d.Wfc =  sqrt(Efc*Ecoeff);
 React_O2_d.N = Npoints;
 React_O2_d.Estep =  (React_O2_d.Emax -  React_O2_d.Emin)/Npoints;

  //outputfile = fopen("xsct/O2_d.dat", "w");

 for ( i = 0; i <= React_O2_d.N; i++ )
 {
 E = i*React_O2_d.Estep + React_O2_d.Emin;

 tmp = 0.;

 if      (E/Ecoeff <= edata[0])  tmp = Sdata[0];
 else if (E/Ecoeff >= edata[12]) tmp = Sdata[12];
     else
   {
     int ind = 0;
    while (E/Ecoeff > edata[ind] ) ind++;
    assert( ind >= 1 && "calculated index cases out of bounds access to Sdata" );
    tmp = Sdata[ind - 1] + (Sdata[ind]- Sdata[ind-1])*(E/Ecoeff - edata[ind-1])/(edata[ind]- edata[ind-1]);
   }

 React_O2_d.CS[i] = CScoeff*tmp*1.e-16*coll_fac_ntrl_ntrl;

 // if (!Node)  fprintf(outputfile, "%f  %15.7e \n", E/Ecoeff, React_O2_d.CS[i]/CScoeff  );
   }
 //fclose(outputfile); 
 
 
// Excitation singlet delta cross section O2 [1e-20 m^2]
//	 Source : 
 

edata[0]   = 0.977; 	Sdata[0]   = 0.;           
edata[1]   = 1.5;      	Sdata[1]   = 0.0058;  
edata[2]   = 2.;        Sdata[2]   = 0.0153;    
edata[3]   = 3.;        Sdata[3]   = 0.038;     
edata[4]   = 3.5;      	Sdata[4]   = 0.049;     
edata[5]   = 4;         Sdata[5]   = 0.057;     
edata[6]   = 5.;        Sdata[6]   = 0.074;     
edata[7]   = 5.62;    	Sdata[7]   = 0.0825;   
edata[8]   = 5.91;   	Sdata[8]   = 0.0862;   
edata[9]   = 6.19;   	Sdata[9]   = 0.0888;    
edata[10]  =6.53;     	Sdata[10]  =0.0908;    
edata[11]  =6.99;     	Sdata[11]  =0.0914;    
edata[12]  =7.61;     	Sdata[12]  =0.0891;    
edata[13]  =7.89;    	Sdata[13]  =0.0863;     
edata[14]  =8.96;    	Sdata[14]  =0.0768;    
edata[15]  =10.04;  	Sdata[15]  =0.0679;    
edata[16]  =13.;       	Sdata[16]  =0.0527;    
edata[17]  =15.1;    	Sdata[17]  =0.0455;    
edata[18]  =17.5;    	Sdata[18]  =0.0387;    
edata[19]  =20.5;    	Sdata[19]  =0.0324;    
edata[20]  =24.9;    	Sdata[20]  =0.0256;    
edata[21]  =30.9;    	Sdata[21]  =0.0196;    
edata[22]  =41.;       	Sdata[22]  =0.0137;    
edata[23]  =45.;        Sdata[23]  =0.012;       

 

 Emin = 0.977;             //eV
 Eth =  0.977;            //eV
 Emax = 45.;         //eV
 Efc =  0.;          //eV


//Ecoeff = (Omega_pe/dx)*(Omega_pe/dx)*(1.+M_ions[13])/M_ions[13]/(0.5*T_e);

 Ecoeff = (dt/dr)*(dt/dr)*(1.+mn_over_me)/mn_over_me/(0.5*T_e0);

 React_O2_ex1D.Emin = Emin*Ecoeff;
 React_O2_ex1D.Eth =  Eth*Ecoeff;
 React_O2_ex1D.Emax = Emax*Ecoeff;
 //React_O2_ex1D.Wfc =  sqrt(Efc*Ecoeff);
 React_O2_ex1D.N = Npoints;
 React_O2_ex1D.Estep =  (React_O2_ex1D.Emax -  React_O2_ex1D.Emin)/Npoints;

 outputfile = fopen("xsct/O2_ex1D.dat", "w");

 for ( i = 0; i <= React_O2_ex1D.N; i++ )
 {
 E = i*React_O2_ex1D.Estep + React_O2_ex1D.Emin;

 tmp = 0.;

 if      (E/Ecoeff <= edata[0])  tmp = Sdata[0];
 else if (E/Ecoeff >= edata[23]) tmp = Sdata[23];
     else
   {
     int ind = 0;
    while (E/Ecoeff > edata[ind] ) ind++;
    assert( ind >= 1 && "calculated index cases out of bounds access to Sdata" );
    tmp = Sdata[ind - 1] + (Sdata[ind]- Sdata[ind-1])*(E/Ecoeff - edata[ind-1])/(edata[ind]- edata[ind-1]);
   }

 React_O2_ex1D.CS[i] = CScoeff*tmp*1.e-16*coll_fac_ntrl_ntrl;

 fprintf(outputfile, "%f  %15.7e \n", E/Ecoeff, React_O2_ex1D.CS[i]/CScoeff  );
   }
 fclose(outputfile); 
 
 

// Excitation singlet delta cross section O2 [1e-?? m^2]  from Brunger et al
//	 Source :  Brunger , ... ,  Physics Reports 357 (2002) 215-458
 

edata[0]   =  0.977 ;   	Sdata[0]   =  0.0        ;     
edata[1]   =  1.1698 ;     	Sdata[1]   =  0.654872   ;    
edata[2]   =  1.33764 ;    	Sdata[2]   =  1.20776     ;      
edata[3]   =  1.60136 ;    	Sdata[3]   =  2.05836     ;     
edata[4]   =  1.84101 ;    	Sdata[4]   =  2.76725     ;     
edata[5]   =  2.02094 ;    	Sdata[5]   =  3.41932     ;    
edata[6]   =  2.21273 ;    	Sdata[6]   =  4.02893     ;     
edata[7]   =  2.36858 ;    	Sdata[7]   =  4.53929     ;     
edata[8]   =  2.54835 ;   	Sdata[8]   =  5.0922       ;    
edata[9]   =  2.69202 ;   	Sdata[9]   =  5.44672     ;      
edata[10]  = 2.97904 ;     	Sdata[10]  = 5.94328     ;      
edata[11]  = 3.31412 ;     	Sdata[11]  = 6.66661     ;     
edata[12]  = 3.61295 ;     	Sdata[12]  = 7.09239     ;      
edata[13]  = 3.97134 ;    	Sdata[13]  = 7.47584     ;       
edata[14]  = 4.3057 ;    	Sdata[14]  = 7.74591     ;     
edata[15]  = 4.69958 ;    	Sdata[15]  = 7.94532     ;      
edata[16]  = 4.99782 ;    	Sdata[16]  = 8.00283     ;      
edata[17]  = 5.30772 ;    	Sdata[17]  = 7.89039     ;     
edata[18]  = 5.534 ;    	Sdata[18]  = 7.69273     ;      
edata[19]  = 5.78406 ;    	Sdata[19]  = 7.45265     ;     
edata[20]  = 5.97448 ;    	Sdata[20]  = 7.19823     ;     
edata[21]  = 6.12927 ;    	Sdata[21]  = 7.04286     ;      
edata[22]  = 6.41511 ;    	Sdata[22]  = 6.80288     ;      
edata[23]  = 6.7368 ;       Sdata[23]  = 6.60549     ;      
edata[24]  = 7.04674 ;    	Sdata[24]  = 6.52138     ;      
edata[25]  = 7.27352 ;    	Sdata[25]  = 6.63534     ;      
edata[26]  = 7.36949 ;     	Sdata[26]  = 6.98972     ;     
edata[27]  = 7.50174 ;    	Sdata[27]  = 7.65581     ;      
edata[28]  = 7.61033 ;    	Sdata[28]  = 8.44932     ;      
edata[29]  = 7.69575 ;     	Sdata[29]  = 9.66769     ;     
edata[30]  = 7.85395 ;    	Sdata[30]  = 11.6511     ;      
edata[31]  = 7.98741 ;     	Sdata[31]  = 13.0821     ;      
edata[32]  = 8.04801 ;    	Sdata[32]  = 13.6914     ;     
edata[33]  = 8.15566 ;     	Sdata[33]  = 13.89         ;      
edata[34]  = 8.34605 ;     	Sdata[34]  = 13.6214     ;      
edata[35]  = 8.58279 ;    	Sdata[35]  = 12.5031     ;     
edata[36]  = 8.87917 ;    	Sdata[36]  = 11.3849     ;       
edata[37]  = 9.30659 ;    	Sdata[37]  = 10.1822     ;       
edata[38]  = 9.7227 ;    	Sdata[38]  = 9.36182     ;     
edata[39]  = 10.2107 ;    	Sdata[39]  = 8.7683       ;      
edata[40]  = 10.7825 ;    	Sdata[40]  = 8.31666     ;       
edata[41]  = 11.6402 ;    	Sdata[41]  = 7.73836     ;       
edata[42]  = 12.3195 ;    	Sdata[42]  = 7.42867     ;     
edata[43]  = 12.9394 ;     	Sdata[43]  = 7.27462     ;      
edata[44]  = 13.5951 ;     	Sdata[44]  = 7.12067     ;      
edata[45]  = 14.3224 ;    	Sdata[45]  = 6.98108     ;     
edata[46]  = 15.0855 ;    	Sdata[46]  = 6.91242     ;       
 

 

 Emin = 0.977;             //eV
 Eth =  0.977;            //eV
 Emax = 15.;         //eV
 Efc =  0.;          //eV


//Ecoeff = (Omega_pe/dx)*(Omega_pe/dx)*(1.+M_ions[13])/M_ions[13]/(0.5*T_e);

 Ecoeff = (dt/dr)*(dt/dr)*(1.+mn_over_me)/mn_over_me/(0.5*T_e0);

 React_O2_ex1Dbrung.Emin = Emin*Ecoeff;
 React_O2_ex1Dbrung.Eth =  Eth*Ecoeff;
 React_O2_ex1Dbrung.Emax = Emax*Ecoeff;
// React_O2_ex1Dbrung.Wfc =  sqrt(Efc*Ecoeff);
 React_O2_ex1Dbrung.N = Npoints;
 React_O2_ex1Dbrung.Estep =  (React_O2_ex1Dbrung.Emax -  React_O2_ex1Dbrung.Emin)/Npoints;

 outputfile = fopen("xsct/O2_ex1D_brung.dat", "w");

 for ( i = 0; i <= React_O2_ex1Dbrung.N; i++ )
 {
 E = i*React_O2_ex1Dbrung.Estep + React_O2_ex1Dbrung.Emin;

 tmp = 0.;

 if      (E/Ecoeff <= edata[0])  tmp = Sdata[0];
 else if (E/Ecoeff >= edata[46]) tmp = Sdata[46];
     else
   {
     int ind = 0;
    while (E/Ecoeff > edata[ind] ) ind++;
    assert( ind >= 1 && "calculated index cases out of bounds access to Sdata" );
    tmp = Sdata[ind - 1] + (Sdata[ind]- Sdata[ind-1])*(E/Ecoeff - edata[ind-1])/(edata[ind]- edata[ind-1]);
   }

 React_O2_ex1Dbrung.CS[i] = CScoeff*tmp*1.e-18*coll_fac_ntrl_ntrl;

  fprintf(outputfile, "%f  %15.7e \n", E/Ecoeff, React_O2_ex1Dbrung.CS[i]/CScoeff  );
   }
 fclose(outputfile); 
 
 
//Excitation singlet sigma cross section O2 [1e-20 m^2]
// Source
 
edata[0]   = 1.627; 	Sdata[0]   = 0;           
edata[1]   = 2;         	Sdata[1]   = 0.0026; 
edata[2]   = 3;         	Sdata[2]   = 0.0097;   
edata[3]   = 3.5;      	Sdata[3]   = 0.0133;  
edata[4]   = 4;         	Sdata[4]   = 0.0149;  
edata[5]   = 5;         	Sdata[5]   = 0.0182;  
edata[6]   = 5.69;    	Sdata[6]   = 0.0194;  
edata[7]   = 6.54;    	Sdata[7]   = 0.0194;  
edata[8]   = 7.34;   	Sdata[8]   = 0.0191;  
edata[9]   = 8.41;   	Sdata[9]   = 0.0183;   
edata[10]  =9.26;     	Sdata[10]  =0.0174;   
edata[11]  =10;        	Sdata[11]  =0.016;     
edata[12]  =13;        	Sdata[12]  =0.013;     
edata[13]  =14.9;    	 Sdata[13]  =0.013;      
edata[14]  =17;       	Sdata[14]  =0.013;     
edata[15]  =19.4;    	Sdata[15]  =0.0125;   
edata[16]  =20.7;    	Sdata[16]  =0.0125;   
edata[17]  =22.5;    	Sdata[17]  =0.011;     
edata[18]  =24;       	Sdata[18]  =0.01;       
edata[19]  =28;       	Sdata[19]  =0.008;     
edata[20]  =35.1;    	Sdata[20]  =0.0063;   
edata[21]  =41.9;    	Sdata[21]  =0.0018;   
edata[22]  =45.1;    	Sdata[22]  =5.e-4;      

 

 Emin = 1.627;             //eV
 Eth =  1.627;            //eV
 Emax = 45.1;         //eV
 Efc =  0.;          //eV


//Ecoeff = (Omega_pe/dx)*(Omega_pe/dx)*(1.+M_ions[13])/M_ions[13]/(0.5*T_e);

 Ecoeff = (dt/dr)*(dt/dr)*(1.+mn_over_me)/mn_over_me/(0.5*T_e0);

 React_O2_ex1S.Emin = Emin*Ecoeff;
 React_O2_ex1S.Eth =  Eth*Ecoeff;
 React_O2_ex1S.Emax = Emax*Ecoeff;
// React_O2_ex1S.Wfc =  sqrt(Efc*Ecoeff);
 React_O2_ex1S.N = Npoints;
 React_O2_ex1S.Estep =  (React_O2_ex1S.Emax -  React_O2_ex1S.Emin)/Npoints;

  outputfile = fopen("xsct/O2_ex1S.dat", "w");


 for ( i = 0; i <= React_O2_ex1S.N; i++ )
 {
 E = i*React_O2_ex1S.Estep + React_O2_ex1S.Emin;

 tmp = 0.;

 if      (E/Ecoeff <= edata[0])  tmp = Sdata[0];
 else if (E/Ecoeff >= edata[22]) tmp = Sdata[22];
     else
   {
     int ind = 0;
    while (E/Ecoeff > edata[ind] ) ind++;
    assert( ind >= 1 && "calculated index cases out of bounds access to Sdata" );
    tmp = Sdata[ind - 1] + (Sdata[ind]- Sdata[ind-1])*(E/Ecoeff - edata[ind-1])/(edata[ind]- edata[ind-1]);
   }

 React_O2_ex1S.CS[i] = CScoeff*tmp*1.e-16*coll_fac_ntrl_ntrl;

 fprintf(outputfile, "%f  %15.7e \n", E/Ecoeff, React_O2_ex1S.CS[i]/CScoeff  );
   }
 fclose(outputfile); 
 
  
 //Excitation singlet sigma cross section O2 [1e-20 m^2] from Brunger data 
// Source:  Brunger , ... ,  Physics Reports 357 (2002) 215-458
/* 
edata[0]   =  1.627 ;   	Sdata[0]   =  0.0 ;     
edata[1]   =  1.86459 ;     Sdata[1]   =  0.179347     ;    
edata[2]   =  2.18488 ;    	Sdata[2]   =  0.44724       ;      
edata[3]   =  2.47973 ;    	Sdata[3]   =  0.736187     ;     
edata[4]   =  2.81314 ;    	Sdata[4]   =  1.08851       ;     
edata[5]   =  3.17206 ;    	Sdata[5]   =  1.43385       ;    
edata[6]   =  3.39001 ;    	Sdata[6]   =  1.65231       ;     
edata[7]   =  3.65902 ;    	Sdata[7]   =  1.87087       ;     
edata[8]   =  4.00458 ;   	Sdata[8]   =  2.07551       ;    
edata[9]   =  4.31144 ;   	Sdata[9]   =  2.18863       ;      
edata[10]  = 4.61813 ;     	Sdata[10]  = 2.25954       ;      
edata[11]  = 5.02688 ;     	Sdata[11]  = 2.31658       ;     
edata[12]  = 5.33326 ;     	Sdata[12]  = 2.31716       ;      
edata[13]  = 5.7416 ;    	Sdata[13]  = 2.27572       ;       
edata[14]  = 6.03498 ;    	Sdata[14]  = 2.22             ;     
edata[15]  = 6.27729 ;    	Sdata[15]  = 2.16418       ;      
edata[16]  = 6.46867 ;    	Sdata[16]  = 2.1364         ;      
edata[17]  = 6.67265 ;    	Sdata[17]  = 2.07348       ;     
edata[18]  = 6.81302 ;    	Sdata[18]  = 2.05967       ;      
edata[19]  = 6.9666 ;    	Sdata[19]  = 2.1514         ;     
edata[20]  = 7.12069 ;    	Sdata[20]  = 2.36271       ;     
edata[21]  = 7.23723 ;    	Sdata[21]  = 2.7498         ;      
edata[22]  = 7.35431 ;    	Sdata[22]  = 3.2635         ;      
edata[23]  = 7.51098 ;      Sdata[23]  = 4.07974       ;      
edata[24]  = 7.64202 ;    	Sdata[24]  = 4.87483       ;      
edata[25]  = 7.78529 ;    	Sdata[25]  = 5.54333       ;      
edata[26]  = 7.88898 ;     	Sdata[26]  = 5.90929       ;     
edata[27]  = 7.96632 ;    	Sdata[27]  = 6.08528       ;      
edata[28]  = 8.04304 ;    	Sdata[28]  = 6.11356       ;      
edata[29]  = 8.195 ;     	Sdata[29]  = 5.82545       ;     
edata[30]  = 8.3957 ;    	Sdata[30]  = 4.98878       ;      
edata[31]  = 8.52135 ;     	Sdata[31]  = 4.51774       ;      
edata[32]  = 8.68546 ;    	Sdata[32]  = 4.08194       ;     
edata[33]  = 8.95157 ;     	Sdata[33]  = 3.6182         ;      
edata[34]  = 9.34594 ;     	Sdata[34]  = 3.29537       ;      
edata[35]  = 9.85514 ;    	Sdata[35]  = 2.9587         ;     
edata[36]  = 10.5814 ;    	Sdata[36]  = 2.6365         ;       
edata[37]  = 11.2315 ;    	Sdata[37]  = 2.41263       ;       
edata[38]  = 12.0734 ;    	Sdata[38]  = 2.24539       ;     
edata[39]  = 12.9921 ;    	Sdata[39]  = 2.14864       ;      
edata[40]  = 14.051 ;    	Sdata[40]  = 1.99588       ;       
edata[41]  = 14.7783 ;    	Sdata[41]  = 1.89876       ;       
edata[42]  = 15.1101 ;    	Sdata[42]  = 1.86422       ;     


 

 Emin = 1.627;             //eV
 Eth =  1.627;            //eV
 Emax = 15.1;         //eV
 Efc =  0.;          //eV


//Ecoeff = (Omega_pe/dx)*(Omega_pe/dx)*(1.+M_ions[13])/M_ions[13]/(0.5*T_e);

 Ecoeff = (dt/dr)*(dt/dr)*(1.+mn_over_me)/mn_over_me/(0.5*T_e0);

 React_O2_ex1Sbrung.Emin = Emin*Ecoeff;
 React_O2_ex1Sbrung.Eth =  Eth*Ecoeff;
 React_O2_ex1Sbrung.Emax = Emax*Ecoeff;
// React_O2_ex1Sbrung.Wfc =  sqrt(Efc*Ecoeff);
 React_O2_ex1Sbrung.N = Npoints;
 React_O2_ex1Sbrung.Estep =  (React_O2_ex1Sbrung.Emax -  React_O2_ex1Sbrung.Emin)/Npoints;

   outputfile = fopen("xsct/O2_ex1S_brung.dat", "w");

 ind = 0;

 for ( i = 0; i <= React_O2_ex1Sbrung.N; i++ )
 {
 E = i*React_O2_ex1Sbrung.Estep + React_O2_ex1Sbrung.Emin;

 tmp = 0.;

 if      (E/Ecoeff <= edata[0])  tmp = Sdata[0];
 else if (E/Ecoeff >= edata[42]) tmp = Sdata[42];
     else
   {
    while (E/Ecoeff > edata[ind] ) ind++;
    tmp = Sdata[ind - 1] + (Sdata[ind]- Sdata[ind-1])*(E/Ecoeff - edata[ind-1])/(edata[ind]- edata[ind-1]);
   }

 React_O2_ex1Sbrung.CS[i] = CScoeff*tmp*1.e-18;

 fprintf(outputfile, "%f  %15.7e \n", E/Ecoeff, React_O2_ex1Sbrung.CS[i]/CScoeff  );
   }
  fclose(outputfile);
 */ 
 
 
 
 
 
//	! Excitation > 4.5 eV cross section O2 [1e-20 m^2]   ?? no idea what is it, knm
//	! Source : 

edata[0]   = 4.5; 		Sdata[0]   = 0.;          
edata[1]   = 4.8;  		Sdata[1]   = 0.003;  
edata[2]   = 5.;    	Sdata[2]   = 0.009;    
edata[3]   = 5.5;  		Sdata[3]   = 0.03;     
edata[4]   = 6.;    	Sdata[4]   = 0.065;   
edata[5]   = 6.5;  		Sdata[5]   = 0.085;   
edata[6]   = 7;     	Sdata[6]   = 0.095;   
edata[7]   = 7.5;  		Sdata[7]   = 0.1;       
edata[8]   = 8.;    	Sdata[8]   = 0.1;       
edata[9]   = 9.;    	Sdata[9]   = 0.085;    
edata[10]  = 10.;    	Sdata[10]  = 0.07;      
edata[11]  = 12.;    	Sdata[11]  = 0.045;    
edata[12]  = 15.;    	Sdata[12]  = 0.;          
 


 Emin = 4.5;             //eV
 Eth =  4.5;            //eV
 Emax = 15.;         //eV
 Efc =  0.;          //eV


//Ecoeff = (Omega_pe/dx)*(Omega_pe/dx)*(1.+M_ions[13])/M_ions[13]/(0.5*T_e);

 Ecoeff = (dt/dr)*(dt/dr)*(1.+mn_over_me)/mn_over_me/(0.5*T_e0);

 React_O2_ex_g45.Emin = Emin*Ecoeff;
 React_O2_ex_g45.Eth =  Eth*Ecoeff;
 React_O2_ex_g45.Emax = Emax*Ecoeff;
// React_O2_ex_g45.Wfc =  sqrt(Efc*Ecoeff);
 React_O2_ex_g45.N = Npoints;
 React_O2_ex_g45.Estep =  (React_O2_ex_g45.Emax -  React_O2_ex_g45.Emin)/Npoints;

  //outputfile = fopen("xsct/O2_ex_g45.dat", "w");

 for ( i = 0; i <= React_O2_ex_g45.N; i++ )
 {
 E = i*React_O2_ex_g45.Estep + React_O2_ex_g45.Emin;

 tmp = 0.;

 if      (E/Ecoeff <= edata[0])  tmp = Sdata[0];
 else if (E/Ecoeff >= edata[12]) tmp = Sdata[12];
     else
   {
     int ind = 0;
    while (E/Ecoeff > edata[ind] ) ind++;
    assert( ind >= 1 && "calculated index cases out of bounds access to Sdata" );
    tmp = Sdata[ind - 1] + (Sdata[ind]- Sdata[ind-1])*(E/Ecoeff - edata[ind-1])/(edata[ind]- edata[ind-1]);
   }

 React_O2_ex_g45.CS[i] = CScoeff*tmp*1.e-16*coll_fac_ntrl_ntrl;

  //if (!Node)  fprintf(outputfile, "%f  %15.7e \n", E/Ecoeff, React_O2_ex_g45.CS[i]/CScoeff  );
   }
 //fclose(outputfile); 
 
 
 
//! Excitation > 6.0 eV cross section O2 [1e-20 m^2]
//	! Source : 
edata[0]   = 6.;     	Sdata[0]   = 0.;          
edata[1]   = 7.;      	Sdata[1]   = 0.15;    
edata[2]   = 7.8;   	Sdata[2]   = 0.23;      
edata[3]   = 9.;      	Sdata[3]   = 0.23;     
edata[4]   = 10.;    	Sdata[4]   = 0.21;     
edata[5]   = 12.;    	Sdata[5]   = 0.165;   
edata[6]   = 15.;    	Sdata[6]   = 0.105;   
edata[7]   = 17.;    	Sdata[7]   = 0.065;   
edata[8]   = 20.;   	Sdata[8]   = 0.0475; 
edata[9]   = 45.;   	Sdata[9]   = 0.019;    

 Emin = 6.;             //eV
 Eth =  6.;            //eV
 Emax = 45.;         //eV
 Efc =  0.;          //eV

//Ecoeff = (Omega_pe/dx)*(Omega_pe/dx)*(1.+M_ions[13])/M_ions[13]/(0.5*T_e);

 Ecoeff = (dt/dr)*(dt/dr)*(1.+mn_over_me)/mn_over_me/(0.5*T_e0);

 React_O2_ex_g6.Emin = Emin*Ecoeff;
 React_O2_ex_g6.Eth =  Eth*Ecoeff;
 React_O2_ex_g6.Emax = Emax*Ecoeff;
// React_O2_ex_g6.Wfc =  sqrt(Efc*Ecoeff);
 React_O2_ex_g6.N = Npoints;
 React_O2_ex_g6.Estep =  (React_O2_ex_g6.Emax -  React_O2_ex_g6.Emin)/Npoints;

  //outputfile = fopen("xsct/O2_ex_g6.dat", "w");

 for ( i = 0; i <= React_O2_ex_g6.N; i++ )
 {
 E = i*React_O2_ex_g6.Estep + React_O2_ex_g6.Emin;

 tmp = 0.;

 if      (E/Ecoeff <= edata[0])  tmp = Sdata[0];
 else if (E/Ecoeff >= edata[9])  tmp = Sdata[9];
     else
   {
     int ind = 0;
    while (E/Ecoeff > edata[ind] ) ind++;
    assert( ind >= 1 && "calculated index cases out of bounds access to Sdata" );
    tmp = Sdata[ind - 1] + (Sdata[ind]- Sdata[ind-1])*(E/Ecoeff - edata[ind-1])/(edata[ind]- edata[ind-1]);
   }

 React_O2_ex_g6.CS[i] = CScoeff*tmp*1.e-16*coll_fac_ntrl_ntrl;

  //if (!Node)  fprintf(outputfile, "%f  %15.7e \n", E/Ecoeff, React_O2_ex_g6.CS[i]/CScoeff  );
   }
 //fclose(outputfile); 
 

 //! Excitation > 8.4 eV cross section O2 [1e-20 m^2]
 //	! Source 
 
edata[0]   = 8.4;    	Sdata[0]   = 0.;       
edata[1]   = 9.4;     	Sdata[1]   = 1.;      
edata[2]   = 30.;      	Sdata[2]   = 0.9;     
edata[3]   = 50.;      	Sdata[3]   = 0.7;    
edata[4]   = 100.;    	Sdata[4]   = 0.54;  
edata[5]   = 150.;    	Sdata[5]   = 0.32;  

 Emin = 8.4;             //eV
 Eth =  8.4;            //eV
 Emax = 150.;         //eV
 Efc =  0.;          //eV

//Ecoeff = (Omega_pe/dx)*(Omega_pe/dx)*(1.+M_ions[13])/M_ions[13]/(0.5*T_e);

 Ecoeff = (dt/dr)*(dt/dr)*(1.+mn_over_me)/mn_over_me/(0.5*T_e0);

 React_O2_ex_g84.Emin = Emin*Ecoeff;
 React_O2_ex_g84.Eth =  Eth*Ecoeff;
 React_O2_ex_g84.Emax = Emax*Ecoeff;
// React_O2_ex_g84.Wfc =  sqrt(Efc*Ecoeff);
 React_O2_ex_g84.N = Npoints;
 React_O2_ex_g84.Estep =  (React_O2_ex_g84.Emax -  React_O2_ex_g84.Emin)/Npoints;

  //outputfile = fopen("xsct/O2_ex_g84.dat", "w");

 for ( i = 0; i <= React_O2_ex_g84.N; i++ )
 {
 E = i*React_O2_ex_g84.Estep + React_O2_ex_g84.Emin;

 tmp = 0.;

 if      (E/Ecoeff <= edata[0])  tmp = Sdata[0];
 else if (E/Ecoeff >= edata[5])  tmp = Sdata[5];
     else
   {
     int ind = 0;
    while (E/Ecoeff > edata[ind] ) ind++;
    assert( ind >= 1 && "calculated index cases out of bounds access to Sdata" );
    tmp = Sdata[ind - 1] + (Sdata[ind]- Sdata[ind-1])*(E/Ecoeff - edata[ind-1])/(edata[ind]- edata[ind-1]);
   }

 React_O2_ex_g84.CS[i] = CScoeff*tmp*1.e-16*coll_fac_ntrl_ntrl;

 // if (!Node)  fprintf(outputfile, "%f  %15.7e \n", E/Ecoeff, React_O2_ex_g84.CS[i]/CScoeff  );
   }
 //fclose(outputfile); 
 

//! Vibration > nu=1 cross section O2 [1e-20 m^2]
//	! Source : 
 
edata[0]   = 0.19;   	Sdata[0]   = 0.;             
edata[1]   = 0.2;      	Sdata[1]   = 1.e-3;       
edata[2]   = 0.21;    	Sdata[2]   = 1.e-3;         
edata[3]   = 0.23;    	Sdata[3]   = 0.;             
edata[4]   = 0.32;    	Sdata[4]   = 0.;             
edata[5]   = 0.33;    	Sdata[5]   = 0.415;      
edata[6]   = 0.35;    	Sdata[6]   = 0.;             
edata[7]   = 0.44;    	Sdata[7]   = 0.;             
edata[8]   = 0.45;   	Sdata[8]   = 1.35;        
edata[9]   = 0.47;   	Sdata[9]   = 0.;              
edata[10]  =0.56;     	Sdata[10]  =0.;              
edata[11]  =0.57;     	Sdata[11]  =1.85;         
edata[12]  =0.59;     	Sdata[12]  =0.;              
edata[13]  =0.68;    	Sdata[13]  =0.;               
edata[14]  =0.69;    	Sdata[14]  =1.65;         
edata[15]  =0.71;    	Sdata[15]  =0.;              
edata[16]  =0.79;    	Sdata[16]  =0.;              
edata[17]  =0.8;      	Sdata[17]  =1.;              
edata[18]  =0.82;    	Sdata[18]  =0.;              
edata[19]  =0.9;      	Sdata[19]  =0.;              
edata[20]  =0.91;    	Sdata[20]  =0.6;           
edata[21]  =0.93;    	Sdata[21]  =0.;              
edata[22]  =1.02;    	Sdata[22]  =0.;              
edata[23]  =1.03;     	Sdata[23]  =0.285;        
edata[24]  =1.05;    	Sdata[24]  =0.;              
edata[25]  =1.13;    	Sdata[25]  =0.;              
edata[26]  =1.14;     	Sdata[26]  =0.1125;     
edata[27]  =1.16;    	Sdata[27]  =0.;              
edata[28]  =1.23;    	Sdata[28]  =0.;              
edata[29]  =1.24;     	Sdata[29]  =0.0475;     
edata[30]  =1.26;    	Sdata[30]  =0.;              
edata[31]  =1.34;     	Sdata[31]  =0.;              
edata[32]  =1.35;    	Sdata[32]  =0.0165;     
edata[33]  =1.37;     	Sdata[33]  =0.;              
edata[34]  =1.44;     	Sdata[34]  =0.;              
edata[35]  =1.45;    	Sdata[35]  =0.0055;     
edata[36]  =1.47;    	Sdata[36]  =0.;               
edata[37]  =1.54;    	Sdata[37]  =0.;               
edata[38]  =1.55;    	Sdata[38]  =0.0019;     
edata[39]  =1.57;    	Sdata[39]  =0.;              
edata[40]  =1.63;    	Sdata[40]  =0.;               
edata[41]  =1.65;    	Sdata[41]  =6.e-4;          
edata[42]  =1.67;    	Sdata[42]  =0.;              
edata[43]  =3.5;       	Sdata[43]  =0.;              
edata[44]  =4.;         Sdata[44]  =0.;              
edata[45]  =5.;         Sdata[45]  =0.042;       
edata[46]  =6.;         Sdata[46]  =0.1;           
edata[47]  =7.;         Sdata[47]  =0.176;        
edata[48]  =8.;         Sdata[48]  =0.231;        
edata[49]  =9.;         Sdata[49]  =0.247;        
edata[50]  =10.;       	Sdata[50]  =0.234;        
edata[51]  =11.;       	Sdata[51]  =0.186;        
edata[52]  =12.;       	Sdata[52]  =0.143;        
edata[53]  =13.;       	Sdata[53]  =0.102;        
edata[54]  =14.;       	Sdata[54]  =0.071;       
edata[55]  =15.;        Sdata[55]  =0.04;         
edata[56]  =20.;        Sdata[56]  =0.01;         
edata[57]  =45.;      	Sdata[57]  =0.;              
 
 
 Emin = 0.19;            //eV
 Eth =  0.19;            //eV
 Emax = 45.;             //eV
 Efc =  0.;              //eV

//Ecoeff = (Omega_pe/dx)*(Omega_pe/dx)*(1.+M_ions[13])/M_ions[13]/(0.5*T_e);

 Ecoeff = (dt/dr)*(dt/dr)*(1.+mn_over_me)/mn_over_me/(0.5*T_e0);

 React_O2_vib_1.Emin = Emin*Ecoeff;
 React_O2_vib_1.Eth =  Eth*Ecoeff;
 React_O2_vib_1.Emax = Emax*Ecoeff;
// React_O2_vib_1.Wfc =  sqrt(Efc*Ecoeff);
 React_O2_vib_1.N = Npoints;
 React_O2_vib_1.Estep =  (React_O2_vib_1.Emax -  React_O2_vib_1.Emin)/Npoints;

 //outputfile = fopen("xsct/O2_vib_1.dat", "w");

 for ( i = 0; i <= React_O2_vib_1.N; i++ )
 {
 E = i*React_O2_vib_1.Estep + React_O2_vib_1.Emin;

 tmp = 0.;

 if      (E/Ecoeff <= edata[0])   tmp = Sdata[0];
 else if (E/Ecoeff >= edata[57])  tmp = Sdata[57];
     else
   {
     int ind = 0;
    while (E/Ecoeff > edata[ind] ) ind++;
    assert( ind >= 1 && "calculated index cases out of bounds access to Sdata" );
    tmp = Sdata[ind - 1] + (Sdata[ind]- Sdata[ind-1])*(E/Ecoeff - edata[ind-1])/(edata[ind]- edata[ind-1]);
   }

 React_O2_vib_1.CS[i] = CScoeff*tmp*1.e-16*coll_fac_ntrl_ntrl;

 // if (!Node)  fprintf(outputfile, "%f  %15.7e \n", E/Ecoeff, React_O2_vib_1.CS[i]/CScoeff  );
   }
 //fclose(outputfile); 
  
 
 
// !=================================
//	! Vibration > nu=2 cross section O2 [1e-20 m^2]
//	! Source : 

edata[0]   = 0.38;   	Sdata[0]   = 0.;           
edata[1]   = 0.44;    	Sdata[1]   = 0.;          
edata[2]   = 0.45;    	Sdata[2]   = 0.;            
edata[3]   = 0.47;    	Sdata[3]   = 0.;           
edata[4]   = 0.56;    	Sdata[4]   = 0.;           
edata[5]   = 0.57;    	Sdata[5]   = 0.14;      
edata[6]   = 0.59;    	Sdata[6]   = 0.;           
edata[7]   = 0.68;    	Sdata[7]   = 0.;           
edata[8]   = 0.69;   	Sdata[8]   = 0.415;    
edata[9]   = 0.71;   	Sdata[9]   = 0.;            
edata[10]  =0.79;     	Sdata[10]  =0.;            
edata[11]  =0.8;       	Sdata[11]  =0.535;     
edata[12]  =0.82;     	Sdata[12]  =0.;            
edata[13]  =0.9;      	Sdata[13]  =0.;             
edata[14]  =0.91;    	Sdata[14]  =0.465;     
edata[15]  =0.93;    	Sdata[15]  =0.;            
edata[16]  =1.02;    	Sdata[16]  =0.;            
edata[17]  =1.03;    	Sdata[17]  =0.315;     
edata[18]  =1.05;    	Sdata[18]  =0.;            
edata[19]  =1.13;    	Sdata[19]  =0.;            
edata[20]  =1.14;    	Sdata[20]  =0.2;         
edata[21]  =1.16;    	Sdata[21]  =0.;            
edata[22]  =1.23;    	Sdata[22]  =0.;            
edata[23]  =1.24;     	Sdata[23]  =0.095;      
edata[24]  =1.26;    	Sdata[24]  =0.;            
edata[25]  =1.34;    	Sdata[25]  =0.;            
edata[26]  =1.35;     	Sdata[26]  =0.04;       
edata[27]  =1.37;    	Sdata[27]  =0.;            
edata[28]  =1.44;    	Sdata[28]  =0.;            
edata[29]  =1.45;     	Sdata[29]  =0.0185;   
edata[30]  =1.47;    	Sdata[30]  =0.;            
edata[31]  =1.54;     	Sdata[31]  =0.;            
edata[32]  =1.55;    	Sdata[32]  =0.0085;   
edata[33]  =1.57;     	Sdata[33]  =0.;            
edata[34]  =1.63;     	Sdata[34]  =0.;            
edata[35]  =1.65;    	Sdata[35]  =0.0034;   
edata[36]  =1.67;    	Sdata[36]  =0.;             
edata[37]  =3.5;      	Sdata[37]  =0.;             
edata[38]  =4.;         Sdata[38]  =0.;            
edata[39]  =5.;         Sdata[39]  =0.028;     
edata[40]  =6.;         Sdata[40]  =0.04;        
edata[41]  =7.;         Sdata[41]  =0.073;      
edata[42]  =8.;         Sdata[42]  =0.094;     
edata[43]  =9.;         Sdata[43]  =0.11;       
edata[44]  =10.;        Sdata[44]  =0.109;     
edata[45]  =11.;        Sdata[45]  =0.093;     
edata[46]  =12.;        Sdata[46]  =0.073;     
edata[47]  =13.;       	Sdata[47]  =0.051;      
edata[48]  =14.;       	Sdata[48]  =0.028;      
edata[49]  =15.;       	Sdata[49]  =0.013;      
edata[50]  =20.;       	Sdata[50]  =0.005;      
edata[51]  =45.;       	Sdata[51]  =0.;             

 Emin = 0.38;            //eV
 Eth =  0.38;            //eV
 Emax = 45.;             //eV
 Efc =  0.;              //eV

//Ecoeff = (Omega_pe/dx)*(Omega_pe/dx)*(1.+M_ions[13])/M_ions[13]/(0.5*T_e);

 Ecoeff = (dt/dr)*(dt/dr)*(1.+mn_over_me)/mn_over_me/(0.5*T_e0);

 React_O2_vib_2.Emin = Emin*Ecoeff;
 React_O2_vib_2.Eth =  Eth*Ecoeff;
 React_O2_vib_2.Emax = Emax*Ecoeff;
// React_O2_vib_2.Wfc =  sqrt(Efc*Ecoeff);
 React_O2_vib_2.N = Npoints;
 React_O2_vib_2.Estep =  (React_O2_vib_2.Emax -  React_O2_vib_2.Emin)/Npoints;

 //outputfile = fopen("xsct/O2_vib_2.dat", "w");


 for ( i = 0; i <= React_O2_vib_2.N; i++ )
 {
 E = i*React_O2_vib_2.Estep + React_O2_vib_2.Emin;

 tmp = 0.;

 if      (E/Ecoeff <= edata[0])   tmp = Sdata[0];
 else if (E/Ecoeff >= edata[51])  tmp = Sdata[51];
     else
   {
     int ind = 0;
    while (E/Ecoeff > edata[ind] ) ind++;
    assert( ind >= 1 && "calculated index cases out of bounds access to Sdata" );
    tmp = Sdata[ind - 1] + (Sdata[ind]- Sdata[ind-1])*(E/Ecoeff - edata[ind-1])/(edata[ind]- edata[ind-1]);
   }

 React_O2_vib_2.CS[i] = CScoeff*tmp*1.e-16*coll_fac_ntrl_ntrl;

 // if (!Node)  fprintf(outputfile, "%f  %15.7e \n", E/Ecoeff, React_O2_vib_2.CS[i]/CScoeff  );
   }
 //fclose(outputfile); 
  
 
//! Vibration > nu=3 cross section O2 [1e-20 m^2]
//	! Source :  


edata[0]   = 0.57;       	Sdata[0]   = 0.;           
edata[1]   = 0.68;        	Sdata[1]   = 0.;          
edata[2]   = 0.69;        	Sdata[2]   = 0.0037;   
edata[3]   = 0.71;        	Sdata[3]   = 0.;           
edata[4]   = 0.79;        	Sdata[4]   = 0.;           
edata[5]   = 0.8;          	Sdata[5]   = 0.0215;  
edata[6]   = 0.82;        	Sdata[6]   = 0.;           
edata[7]   = 0.9;          	Sdata[7]   = 0.;           
edata[8]   = 0.91;       	Sdata[8]   = 0.09;      
edata[9]   = 0.93;       	Sdata[9]   = 0.;            
edata[10]  =1.02;         	Sdata[10]  =0.;            
edata[11]  =1.03;         	Sdata[11]  =0.12;       
edata[12]  =1.05;         	Sdata[12]  =0.;            
edata[13]  =1.13;        	Sdata[13]  =0.;             
edata[14]  =1.14;        	Sdata[14]  =0.115;     
edata[15]  =1.16;        	Sdata[15]  =0.;            
edata[16]  =1.23;        	Sdata[16]  =0.;            
edata[17]  =1.24;        	Sdata[17]  =0.095;     
edata[18]  =1.26;        	Sdata[18]  =0.;            
edata[19]  =1.34;        	Sdata[19]  =0.;            
edata[20]  =1.35;        	Sdata[20]  =0.055;     
edata[21]  =1.37;        	Sdata[21]  =0.;            
edata[22]  =1.44;        	Sdata[22]  =0.;            
edata[23]  =1.45;         	Sdata[23]  =0.03;        
edata[24]  =1.47;        	Sdata[24]  =0.;            
edata[25]  =1.54;        	Sdata[25]  =0.;            
edata[26]  =1.55;         	Sdata[26]  =0.0165;   
edata[27]  =1.57;        	Sdata[27]  =0.;            
edata[28]  =1.63;        	Sdata[28]  =0.;            
edata[29]  =1.65;         	Sdata[29]  =0.008;     
edata[30]  =1.67;        	Sdata[30]  =0.;            
edata[31]  =3.5;           	Sdata[31]  =0.;            
edata[32]  =4.;             Sdata[32]  =0.;            
edata[33]  =5.;             Sdata[33]  =0.;            
edata[34]  =6.;             Sdata[34]  =0.0125;   
edata[35]  =7.;             Sdata[35]  =0.0363;   
edata[36]  =8.;             Sdata[36]  =0.0588;    
edata[37]  =9.;             Sdata[37]  =0.075;      
edata[38]  =10.;           	Sdata[38]  =0.0675;   
edata[39]  =11.;           	Sdata[39]  =0.0563;   
edata[40]  =12.;           	Sdata[40]  =0.0475;    
edata[41]  =13.;           	Sdata[41]  =0.03;        
edata[42]  =14.;           	Sdata[42]  =0.0175;   
edata[43]  =15.;            Sdata[43]  =0.0088;   
edata[44]  =20.;            Sdata[44]  =0.;            
       


 Emin = 0.57;            //eV
 Eth =  0.57;            //eV
 Emax = 20.;             //eV
 Efc =  0.;              //eV

//Ecoeff = (Omega_pe/dx)*(Omega_pe/dx)*(1.+M_ions[13])/M_ions[13]/(0.5*T_e);

 Ecoeff = (dt/dr)*(dt/dr)*(1.+mn_over_me)/mn_over_me/(0.5*T_e0);

 React_O2_vib_3.Emin = Emin*Ecoeff;
 React_O2_vib_3.Eth =  Eth*Ecoeff;
 React_O2_vib_3.Emax = Emax*Ecoeff;
// React_O2_vib_3.Wfc =  sqrt(Efc*Ecoeff);
 React_O2_vib_3.N = Npoints;
 React_O2_vib_3.Estep =  (React_O2_vib_3.Emax -  React_O2_vib_3.Emin)/Npoints;

  //outputfile = fopen("xsct/O2_vib_3.dat", "w");

 for ( i = 0; i <= React_O2_vib_3.N; i++ )
 {
 E = i*React_O2_vib_3.Estep + React_O2_vib_3.Emin;

 tmp = 0.;

 if      (E/Ecoeff <= edata[0])   tmp = Sdata[0];
 else if (E/Ecoeff >= edata[44])  tmp = Sdata[44];
     else
   {
     int ind = 0;
    while (E/Ecoeff > edata[ind] ) ind++;
    assert( ind >= 1 && "calculated index cases out of bounds access to Sdata" );
    tmp = Sdata[ind - 1] + (Sdata[ind]- Sdata[ind-1])*(E/Ecoeff - edata[ind-1])/(edata[ind]- edata[ind-1]);
   }

 React_O2_vib_3.CS[i] = CScoeff*tmp*1.e-16*coll_fac_ntrl_ntrl;

  //if (!Node)  fprintf(outputfile, "%f  %15.7e \n", E/Ecoeff, React_O2_vib_3.CS[i]/CScoeff  );
   }
  //fclose(outputfile); 
     
 
 // ! Vibration > nu=4 cross section O2 [1e-20 m^2]
//	! Source : 

edata[0]   = 0.75;      	Sdata[0]   = 0.;           
edata[1]   = 0.79;       	Sdata[1]   = 0.;          
edata[2]   = 0.8;         	Sdata[2]   = 0.0015;   
edata[3]   = 0.82;       	Sdata[3]   = 0.;           
edata[4]   = 0.9;         	Sdata[4]   = 0.;           
edata[5]   = 0.91;       	Sdata[5]   = 0.0055;  
edata[6]   = 0.93;       	Sdata[6]   = 0.;           
edata[7]   = 1.02;       	Sdata[7]   = 0.;           
edata[8]   = 1.03;      	Sdata[8]   = 3.e-4;     
edata[9]   = 1.05;      	Sdata[9]   = 0.;            
edata[10]  =1.13;        	Sdata[10]  =0.;            
edata[11]  =1.14;        	Sdata[11]  =0.0165;   
edata[12]  =1.16;        	Sdata[12]  =0.;            
edata[13]  =1.23;       	Sdata[13]  =0.;             
edata[14]  =1.24;       	Sdata[14]  =0.0315;   
edata[15]  =1.26;       	Sdata[15]  =0.;            
edata[16]  =1.34;       	Sdata[16]  =0.;            
edata[17]  =1.35;       	Sdata[17]  =0.0335;   
edata[18]  =1.37;       	Sdata[18]  =0.;            
edata[19]  =1.44;       	Sdata[19]  =0.;            
edata[20]  =1.45;       	Sdata[20]  =0.0285;   
edata[21]  =1.47;       	Sdata[21]  =0.;            
edata[22]  =1.54;       	Sdata[22]  =0.;            
edata[23]  =1.55;        	Sdata[23]  =0.0215;    
edata[24]  =1.57;       	Sdata[24]  =0.;            
edata[25]  =1.63;       	Sdata[25]  =0.;            
edata[26]  =1.65;        	Sdata[26]  =0.0165;   
edata[27]  =1.67;       	Sdata[27]  =0.;            
edata[28]  =6.;            	Sdata[28]  =0.;            
edata[29]  =7.;             Sdata[29]  =0.0275;   
edata[30]  =8.;            	Sdata[30]  =0.035;     
edata[31]  =9.;             Sdata[31]  =0.0413;   
edata[32]  =10.;          	Sdata[32]  =0.0462;   
edata[33]  =11.;           	Sdata[33]  =0.0313;   
edata[34]  =12.;           	Sdata[34]  =0.025;     
edata[35]  =13.;          	Sdata[35]  =0.0175;   
edata[36]  =14.;          	Sdata[36]  =0.0088;    
edata[37]  =15.;          	Sdata[37]  =0.;             
          


 Emin = 0.75;            //eV
 Eth =  0.75;            //eV
 Emax = 15.;             //eV
 Efc =  0.;              //eV

//Ecoeff = (Omega_pe/dx)*(Omega_pe/dx)*(1.+M_ions[13])/M_ions[13]/(0.5*T_e);

 Ecoeff = (dt/dr)*(dt/dr)*(1.+mn_over_me)/mn_over_me/(0.5*T_e0);

 React_O2_vib_4.Emin = Emin*Ecoeff;
 React_O2_vib_4.Eth =  Eth*Ecoeff;
 React_O2_vib_4.Emax = Emax*Ecoeff;
// React_O2_vib_4.Wfc =  sqrt(Efc*Ecoeff);
 React_O2_vib_4.N = Npoints;
 React_O2_vib_4.Estep =  (React_O2_vib_4.Emax -  React_O2_vib_4.Emin)/Npoints;

  //outputfile = fopen("xsct/O2_vib_4.dat", "w");

 for ( i = 0; i <= React_O2_vib_4.N; i++ )
 {
 E = i*React_O2_vib_4.Estep + React_O2_vib_4.Emin;

 tmp = 0.;

 if      (E/Ecoeff <= edata[0])   tmp = Sdata[0];
 else if (E/Ecoeff >= edata[37])  tmp = Sdata[37];
     else
   {
     int ind = 0;
    while (E/Ecoeff > edata[ind] ) ind++;
    assert( ind >= 1 && "calculated index cases out of bounds access to Sdata" );
    tmp = Sdata[ind - 1] + (Sdata[ind]- Sdata[ind-1])*(E/Ecoeff - edata[ind-1])/(edata[ind]- edata[ind-1]);
   }

 React_O2_vib_4.CS[i] = CScoeff*tmp*1.e-16*coll_fac_ntrl_ntrl;

 // if (!Node)  fprintf(outputfile, "%f  %15.7e \n", E/Ecoeff, React_O2_vib_4.CS[i]/CScoeff  );
   }
 //fclose(outputfile); 
     
// O2+  on O2 elastic from David
//    up to E=0.502 eV: D.R. Gray and J.A. Rees, J. Phys. B: Atom. Mplec. Phys., Vol. 5, May (1972) 1048.
//     from E= 17 eV: R.F. Stebbings, Ben. R. Turner and A.C.H. Smith, J. Chem. Phys., 38, 9, (1963) 2277.   

edata[0]   =     0.000813/2.;       Sdata[0]    =   166.; //667.5;    // artificial 
edata[1]   =     0.00537/2.;        Sdata[1]    =   166.; //256.9; 
edata[2]   =     0.013/2.;          Sdata[2]    =   166.; 
edata[3]   =     0.0336/2.;         Sdata[3]    =   112.3; 
edata[4]   =     0.0502/2.;         Sdata[4]    =    94.0; 
edata[5]   =     0.0933/2.;         Sdata[5]    =    75.8;
edata[6]   =     0.150/2.;          Sdata[6]    =    63.0; 
edata[7]   =     0.182/2.;          Sdata[7]    =    59.3;
edata[8]   =     0.243/2.;          Sdata[8]    =    54.4; 
edata[9]   =     0.302/2.;          Sdata[9]    =    51.6; 
edata[10]   =    0.383/2.;          Sdata[10]    =   46.8; 
edata[11]   =    0.457/2.;          Sdata[11]    =   44.4; 
edata[12]   =    0.502/2.;          Sdata[12]    =   45.1; 
edata[13]   =   17./2.;             Sdata[13]    =   26.0; 
edata[14]   =   23./2.;             Sdata[14]    =   24.0; 
edata[15]   =   30./2.;             Sdata[15]    =   23.0;
edata[16]   =   40./2.;             Sdata[16]    =   22.0; 
edata[17]   =   50./2.;             Sdata[17]    =   21.0;
edata[18]   =   63./2.;             Sdata[18]    =   20.0; 
edata[19]   =  100./2.;             Sdata[19]    =   19.0;
edata[20]   =  400./2.;             Sdata[20]    =   16.0; 
edata[21]   = 1000./2.;             Sdata[21]    =   15.0; 


 Emin = 0.0;      //eV !! was 0 for paper  0.025 - for 80pa
 Eth = 0;       //eV
 Emax = 200.;    //eV
 Efc = 0.;

// Ecoeff = (Omega_pe/dx)*(Omega_pe/dx)*(M_ions[6]+M_ions[13])/M_ions[13]/M_ions[6]/(0.5*T_e);
 Ecoeff = (dt/dr)*(dt/dr)*(1.+mn_over_me)/mn_over_me/(0.5*T_e0);
 Ecoeff *= dt_ntrl*dt_ntrl;

 CScoeff = collision_amplification_factor * ScaleF* n_e0*dr_0*ncoll_n/Ncell1 / dt_ntrl;      
 
 React_O2p_O2_el.Emin = Emin*Ecoeff;
 React_O2p_O2_el.Eth =  Eth*Ecoeff;
 React_O2p_O2_el.Emax = Emax*Ecoeff;
// React_O2p_O2_el.Wfc =  sqrt(Efc*Ecoeff);

 React_O2p_O2_el.N = Npoints;
 React_O2p_O2_el.Estep =  (React_O2p_O2_el.Emax -  React_O2p_O2_el.Emin)/React_O2p_O2_el.N;

    // Now - calculate coefficients 2 go 2 Dim.less units !!
    //outputfile = fopen("xsct/arp_ar_el.dat", "w");

outputfile = fopen("xsct/O2p_O2_el.dat", "w");

 for ( i = 0; i <= React_O2p_O2_el.N; i++ )
 {
 E = i*React_O2p_O2_el.Estep + React_O2p_O2_el.Emin;

 tmp = 0.;

 if      (E/Ecoeff <= edata[0])   tmp = Sdata[0];
 else if (E/Ecoeff >= edata[20])  tmp = Sdata[20];
     else
   {
     int ind = 0;
    while (E/Ecoeff > edata[ind] ) ind++;
    assert( ind >= 1 && "calculated index cases out of bounds access to Sdata" );
    tmp = Sdata[ind - 1] + (Sdata[ind]- Sdata[ind-1])*(E/Ecoeff - edata[ind-1])/(edata[ind]- edata[ind-1]);
   }

 React_O2p_O2_el.CS[i] = CScoeff*tmp*1e-16*0.5*coll_fac_ntrl_ntrl*ampl_pressure; // Note the factor 0.5!! Because mom transfer CS //remove factor 1*10^(-16) and add /dt_ion paulm 22.6., compare with Xe Ar
  
fprintf(outputfile, "%f  %15.7e \n", E/Ecoeff, React_O2p_O2_el.CS[i]/CScoeff);
 }
fclose(outputfile);


/***********************************************************

        Dissociative attachment //  e + O2 -> O- + O

		Phelps and Kajita compillation
************************************************************/

edata[0]   =  4.2000000e+000  ;        Sdata[0]    =   0.0    ; 
edata[1]   =  4.3000000e+000  ;        Sdata[1]    =   8.7914643e-025    ; 
edata[2]   =  4.4000000e+000  ;        Sdata[2]    =   2.6374393e-024    ; 
edata[3]   =  4.5000000e+000  ;        Sdata[3]    =   4.3957321e-024    ; 
edata[4]   =  4.6000000e+000  ;        Sdata[4]    =   7.0331714e-024    ; 
edata[5]   =  4.7000000e+000  ;        Sdata[5]    =   9.6706107e-024    ; 
edata[6]   =  4.8000000e+000  ;        Sdata[6]    =   1.3187196e-023    ; 
edata[7]   =  4.9000000e+000  ;        Sdata[7]    =   1.7582929e-023    ; 
edata[8]   =  5.0000000e+000  ;        Sdata[8]    =   2.1978661e-023    ; 
edata[9]   =  5.1000000e+000  ;        Sdata[9]    =   2.9011832e-023    ; 
edata[10]  =  5.2000000e+000  ;        Sdata[10]   =   3.6045004e-023    ; 
edata[11]  =  5.3000000e+000  ;        Sdata[11]   =   4.4836468e-023    ; 
edata[12]  =  5.4000000e+000  ;        Sdata[12]   =   5.3627932e-023    ; 
edata[13]  =  5.5000000e+000  ;        Sdata[13]   =   6.3298543e-023    ; 
edata[14]  =  5.6000000e+000  ;        Sdata[14]   =   7.4727447e-023    ; 
edata[15]  =  5.7000000e+000  ;        Sdata[15]   =   8.5277204e-023    ; 
edata[16]  =  5.8000000e+000  ;        Sdata[16]   =   9.5826961e-023    ; 
edata[17]  =  5.9000000e+000  ;        Sdata[17]   =   1.0461843e-022    ; 
edata[18]  =  6.0000000e+000  ;        Sdata[18]   =   1.1428904e-022    ; 
edata[19]  =  6.1000000e+000  ;        Sdata[19]   =   1.2308050e-022    ; 
edata[20]  =  6.2000000e+000  ;        Sdata[20]   =   1.3099282e-022    ; 
edata[21]  =  6.3000000e+000  ;        Sdata[21]   =   1.3626770e-022    ; 
edata[22]  =  6.4000000e+000  ;        Sdata[22]   =   1.3099282e-022    ; 
edata[23]  =  6.5000000e+000  ;        Sdata[23]   =   1.4066343e-022    ; 
edata[24]  =  6.6000000e+000  ;        Sdata[24]   =   1.3978428e-022    ; 
edata[25]  =  6.7000000e+000  ;        Sdata[25]   =   1.3714684e-022    ; 
edata[26]  =  6.8000000e+000  ;        Sdata[26]   =   1.3363026e-022  ; 
edata[27]  =  6.9000000e+000  ;        Sdata[27]   =   1.2835538e-022  ; 
edata[28]  =  7.0000000e+000  ;        Sdata[28]   =   1.2220135e-022  ; 
edata[29]  =  7.1000000e+000  ;        Sdata[29]   =   1.1428904e-022  ; 
edata[30]  =  7.2000000e+000  ;        Sdata[30]   =   1.0637672e-022  ; 
edata[31]  =  7.3000000e+000  ;        Sdata[31]   =   9.8464400e-023  ; 
edata[32]  =  7.4000000e+000  ;        Sdata[32]   =   8.9672936e-023  ; 
edata[33]  =  7.5000000e+000  ;        Sdata[33]   =   8.1760618e-023  ; 
edata[34]  =  7.6000000e+000  ;        Sdata[34]   =   7.3848300e-023  ; 
edata[35]  =  7.7000000e+000  ;        Sdata[35]   =   6.4177689e-023  ; 
edata[36]  =  7.8000000e+000  ;        Sdata[36]   =   5.7144518e-023  ; 
edata[37]  =  7.9000000e+000  ;        Sdata[37]   =   5.0111346e-023  ; 
edata[38]  =  8.0000000e+000  ;        Sdata[38]   =   4.4836468e-023  ; 
edata[39]  =  8.1000000e+000  ;        Sdata[39]   =   3.8682443e-023  ; 
edata[40]  =  8.2000000e+000  ;        Sdata[40]   =   3.3407564e-023  ; 
edata[41]  =  8.3000000e+000  ;        Sdata[41]   =   2.8132686e-023  ; 
edata[42]  =  8.4000000e+000  ;        Sdata[42]   =   2.8132686e-023  ; 
edata[43]  =  8.5000000e+000  ;        Sdata[43]   =   2.0220368e-023  ; 
edata[44]  =  8.6000000e+000  ;        Sdata[44]   =   1.6703782e-023  ; 
edata[45]  =  8.7000000e+000  ;        Sdata[45]   =   1.4066343e-023  ; 
edata[46]  =  8.8000000e+000  ;        Sdata[46]   =   1.2308050e-023  ; 
edata[47]  =  8.9000000e+000  ;        Sdata[47]   =   1.0549757e-023  ; 
edata[48]  =  9.0000000e+000  ;        Sdata[48]   =   8.7914643e-024  ; 
edata[49]  =  9.1000000e+000  ;        Sdata[49]   =   7.0331714e-024  ; 
edata[50]  =  9.2000000e+000  ;        Sdata[50]   =   7.0331714e-024  ; 
edata[51]  =  9.3000000e+000  ;        Sdata[51]   =   6.1540250e-024  ; 
edata[52]  =  9.4000000e+000  ;        Sdata[52]   =   5.2748786e-024  ; 
edata[53]  =  9.5000000e+000  ;        Sdata[53]   =   4.3957321e-024  ; 
edata[54]  =  9.6000000e+000  ;        Sdata[54]   =   4.3957321e-024  ; 
edata[55]  =  9.8000000e+000  ;        Sdata[55]   =   3.5165857e-024  ; 
edata[56]  =  9.9000000e+000  ;        Sdata[56]   =   3.5165857e-024  ; 
edata[57]  =  1.5000000e+001  ;        Sdata[57]   =   0.0000000e+000  ; 
edata[58]  =  1.6000000e+001  ;        Sdata[58]   =   0.0000000e+000  ; 
edata[59]  =  1.7000000e+001  ;        Sdata[59]   =   0.0000000e+000  ; 
edata[60]  =  1.8000000e+001  ;        Sdata[60]   =   1.7582929e-024  ; 
edata[61]  =  1.9000000e+001  ;        Sdata[61]   =   4.3957321e-024  ; 
edata[62]  =  2.0000000e+001  ;        Sdata[62]   =   6.1540250e-024  ; 
edata[63]  =  2.1000000e+001  ;        Sdata[63]   =   9.6706107e-024  ; 
edata[64]  =  2.2000000e+001  ;        Sdata[64]   =   1.7582929e-023  ; 
edata[65]  =  2.3000000e+001  ;        Sdata[65]   =   2.1978661e-023  ; 
edata[66]  =  2.4000000e+001  ;        Sdata[66]   =   2.9011832e-023  ; 
edata[67]  =  2.5000000e+001  ;        Sdata[67]   =   3.3407564e-023  ; 
edata[68]  =  2.6000000e+001  ;        Sdata[68]   =   3.5165857e-023  ; 
edata[69]  =  2.7000000e+001  ;        Sdata[69]   =   3.6045004e-023  ; 
edata[70]  =  2.8000000e+001  ;        Sdata[70]   =   3.8682443e-023  ; 
edata[71]  =  3.2000000e+001  ;        Sdata[71]   =   4.4836468e-023  ; 
edata[72]  =  3.3000000e+001  ;        Sdata[72]   =   4.7473907e-023  ; 
edata[73]  =  3.6000000e+001  ;        Sdata[73]   =   4.5715614e-023  ; 
edata[74]  =  3.8000000e+001  ;        Sdata[74]   =   4.4836468e-023  ; 
edata[75]  =  4.0000000e+001  ;        Sdata[75]   =   4.3957321e-023  ; 
edata[76]  =  4.6000000e+001  ;        Sdata[76]   =   4.2199029e-023  ; 
edata[77]  =  4.8000000e+001  ;        Sdata[77]   =   4.1319882e-023  ; 
edata[78]  =  5.2000000e+001  ;        Sdata[78]   =   4.0440736e-023  ; 
edata[79]  =  5.4000000e+001  ;        Sdata[79]   =   4.0440736e-023  ; 



 Emin = 4.2;             //eV
 Eth =  4.2;            //eV
 Emax = 54.;         //eV
 Efc =  0.;          //eV


//Ecoeff = (Omega_pe/dx)*(Omega_pe/dx)*(1.+M_ions[13])/M_ions[13]/(0.5*T_e);

 Ecoeff = (dt/dr)*(dt/dr)*(1.+mn_over_me)/mn_over_me/(0.5*T_e0);

 React_O2_dat.Emin = Emin*Ecoeff;
 React_O2_dat.Eth =  Eth*Ecoeff;
 React_O2_dat.Emax = Emax*Ecoeff;
// React_O2_dat.Wfc =  sqrt(Efc*Ecoeff);
 React_O2_dat.N = Npoints;
 React_O2_dat.Estep =  (React_O2_dat.Emax -  React_O2_dat.Emin)/Npoints;

  outputfile = fopen("xsct/O2_dat.dat", "w");

 for ( i = 0; i <= React_O2_dat.N; i++ )
 {
 E = i*React_O2_dat.Estep + React_O2_dat.Emin;

 tmp = 0.;

 if      (E/Ecoeff <= edata[0])  tmp = Sdata[0];
 else if (E/Ecoeff >= edata[79]) tmp = Sdata[79];
     else
   {
     int ind = 0;
    while (E/Ecoeff > edata[ind] ) ind++;
    assert( ind >= 1 && "calculated index cases out of bounds access to Sdata" );
    tmp = Sdata[ind - 1] + (Sdata[ind]- Sdata[ind-1])*(E/Ecoeff - edata[ind-1])/(edata[ind]- edata[ind-1]);
   }

 React_O2_dat.CS[i] = CScoeff*tmp*10000*coll_fac_ntrl_ntrl*ampl_pressure; 
 
    fprintf(outputfile, "%f  %15.7e \n", E/Ecoeff, React_O2_dat.CS[i]/CScoeff/ampl_pressure/coll_fac_ntrl_ntrl  );
   }
 fclose(outputfile); 
  


/***********************************************************

        electron impact detachment //  O- + e -> 2e + O 

		from Devid compilation
************************************************************/



edata[0]   =  1.46;        Sdata[0]    =  0.0; 
edata[1]   =  8.0;         Sdata[1]    =  4.0; 
edata[2]   =  12.0;        Sdata[2]    =  4.9; 
edata[3]   =  14.0;        Sdata[3]    =  5.1; 
edata[4]   =  16.0;        Sdata[4]    =  5.3; 
edata[5]   =  30.0;        Sdata[5]    =  5.3; 
edata[6]   =  68.;         Sdata[6]    =  5.3; 
edata[7]   =  98.;         Sdata[7]    =  4.9; 
edata[8]   =  108.;        Sdata[8]    =  4.7; 
edata[9]   =  148.;        Sdata[9]    =  4.09; 
edata[10]  =  194.;        Sdata[10]   =  3.71; 
edata[11]  =  296.;        Sdata[11]   =  2.88; 
edata[12]  =  396.;        Sdata[12]   =  2.28; 
edata[13]  =  495.;        Sdata[13]   =  1.89; 
edata[14]  =  595.;        Sdata[14]   =  1.62; 
edata[15]  =  694.;        Sdata[15]   =  1.45; 
edata[16]  =  794.;        Sdata[16]   =  1.33; 
edata[17]  =  990.;        Sdata[17]   =  1.18;



 Emin = 1.46;             //eV
 Eth =  1.46;            //eV
 Emax = 990.;         //eV
 Efc =  0.;          //eV


//Ecoeff = (Omega_pe/dx)*(Omega_pe/dx)*(1.+M_ions[4])/M_ions[4]/(0.5*T_e);

 Ecoeff = (dt/dr)*(dt/dr)*(1.+mn_over_me)/mn_over_me/(0.5*T_e0);

 React_O_edet.Emin = Emin*Ecoeff;
 React_O_edet.Eth =  Eth*Ecoeff;
 React_O_edet.Emax = Emax*Ecoeff;
// React_O_edet.Wfc =  sqrt(Efc*Ecoeff);
 React_O_edet.N = Npoints;
 React_O_edet.Estep =  (React_O_edet.Emax -  React_O_edet.Emin)/Npoints;

  outputfile = fopen("xsct/O_edet.dat", "w");

 for ( i = 0; i <= React_O_edet.N; i++ )
 {
 E = i*React_O_edet.Estep + React_O_edet.Emin;

 tmp = 0.;

 if      (E/Ecoeff <= edata[0])  tmp = Sdata[0];
 else if (E/Ecoeff >= edata[17]) tmp = Sdata[17];
     else
   {
     int ind = 0;
    while (E/Ecoeff > edata[ind] ) ind++;
    assert( ind >= 1 && "calculated index cases out of bounds access to Sdata" );
    tmp = Sdata[ind - 1] + (Sdata[ind]- Sdata[ind-1])*(E/Ecoeff - edata[ind-1])/(edata[ind]- edata[ind-1]);
   }

 React_O_edet.CS[i] = CScoeff*tmp*1.e-16;

    fprintf(outputfile, "%f  %15.7e \n", E/Ecoeff, React_O_edet.CS[i]/CScoeff  );
   }
 fclose(outputfile); 

/***********************************************************//**
    Elastic colision   // O- + O2 = O- + O2	
    data from Vahedi and David Tshakaya    
    Now replaced with Franz' analytic firmula:
     3.*3.1415*(0.5*3.77*10^(-24))^(2./3.)/(x^(2./3.)) +1.219*10^-15
 **************************************************************/ 

 Emin = 0.005;      //eV
 Eth =  0.0;       //eV
 Emax = 40.005;    //eV
 Efc = 0.;

 Ecoeff = (dt/dr)*(dt/dr)*(1.+mn_over_me)/mn_over_me/(0.5*T_e0);
 Ecoeff *= dt_ion*dt_ion;

 React_Omin_O2_el.Emin = Emin*Ecoeff;
 React_Omin_O2_el.Eth =  Eth*Ecoeff;
 React_Omin_O2_el.Emax = Emax*Ecoeff;
// React_Omin_O2_el.Wfc =  sqrt(Efc*Ecoeff);

 React_Omin_O2_el.N = Npoints;
 React_Omin_O2_el.Estep =  (React_Omin_O2_el.Emax -  React_Omin_O2_el.Emin)/React_Omin_O2_el.N;

    outputfile = fopen("xsct/Omin_O2_el.dat", "w");


 for ( i = 0; i <= React_Omin_O2_el.N; i++ )
 {
   E = i*React_Omin_O2_el.Estep + React_Omin_O2_el.Emin;
   React_Omin_O2_el.CS[i] = CScoeff*(3.*3.1415*pow(0.5*3.77*1.0e-24*Ecoeff/E, 2./3.) + 1.219*1.0e-15)*coll_fac_ntrl_ntrl*ampl_pressure/dt_ion;
   	 
  
  fprintf(outputfile, "%f %15.7e \n", E/Ecoeff, React_Omin_O2_el.CS[i]/CScoeff/coll_fac_ntrl_ntrl/ampl_pressure*dt_ion);

 }
 fclose(outputfile);



 /********************************************************//**

        electron dissociative recombination  //  e + O2+ -> O + O 
		from picture in Vahedi :-)
************************************************************/



edata[0]   =  0.01;          Sdata[0]    =  18.0; 
edata[1]   =  0.133;         Sdata[1]    =  15.0; 
edata[2]   =  0.27;          Sdata[2]    =  13.3; 
edata[3]   =  0.57;          Sdata[3]    =  10.0; 
edata[4]   =  0.8;           Sdata[4]    =  8.0; 
edata[5]   =  1.0;           Sdata[5]    =  6.5; 
edata[6]   =  2.0;           Sdata[6]    =  2.4; 
edata[7]   =  3.0;           Sdata[7]    =  0.9; 
edata[8]   =  5.0;           Sdata[8]    =  0.11; 
edata[9]   =  7.0;           Sdata[9]    =  0.01;
edata[10]  =  9.0;           Sdata[10]   =  0.0; 


 Emin = 0.0;             //eV
 Eth =  0.0;            //eV
 Emax = 9.0;         //eV
 Efc =  12.06;          //eV


//Ecoeff = (Omega_pe/dx)*(Omega_pe/dx)*(1.+M_ions[6])/M_ions[6]/(0.5*T_e);

 Ecoeff = (dt/dr)*(dt/dr)*(1.+mi_over_me)/mi_over_me/(0.5*T_e0);

 React_O2_dissrec.Emin = Emin*Ecoeff;
 React_O2_dissrec.Eth =  Eth*Ecoeff;
 React_O2_dissrec.Emax = Emax*Ecoeff;
// React_O2_dissrec.Wfc =  Efc*Ecoeff;            // square of e - O2+ relative velocity, corresponding to FC energy
 React_O2_dissrec.N = Npoints;
 React_O2_dissrec.Estep =  (React_O2_dissrec.Emax -  React_O2_dissrec.Emin)/Npoints;

  outputfile = fopen("xsct/O_dissrec.dat", "w");


 for ( i = 0; i <= React_O2_dissrec.N; i++ )
 {
 E = i*React_O2_dissrec.Estep + React_O2_dissrec.Emin;

 tmp = 0.;

 if      (E/Ecoeff <= edata[0])  tmp = Sdata[0];
 else if (E/Ecoeff >= edata[10]) tmp = Sdata[10];
     else
   {
     int ind = 0 ;
    while (E/Ecoeff > edata[ind] ) ind++;
    assert( ind >= 1 && "calculated index cases out of bounds access to Sdata" );
    tmp = Sdata[ind - 1] + (Sdata[ind]- Sdata[ind-1])*(E/Ecoeff - edata[ind-1])/(edata[ind]- edata[ind-1]);
   }

 React_O2_dissrec.CS[i] = CScoeff*tmp*1.e-16;  //here no coll_fac_ntrl_ntrl because no neutrals

    fprintf(outputfile, "%f  %15.7e \n", E/Ecoeff, React_O2_dissrec.CS[i]/CScoeff  );
   }
 fclose(outputfile); 

/***********************************************************

        Detachment on neutrals  //O- + O2 = O + O2 + e
		from David compilation
		added low energy part from Franz (analythic)
************************************************************/


edata[0]   =  1.3   ;        Sdata[0]    =   3.e-22    ; 
edata[1]   =  1.8   ;        Sdata[1]    =   9.e-22    ; 
edata[2]   =  2.1   ;        Sdata[2]    =   1.5e-21   ; 
edata[3]   =  2.7   ;        Sdata[3]    =   3.7e-21   ; 
edata[4]   =  3.4   ;        Sdata[4]    =   5.8e-21   ; 
edata[5]   =  4.0   ;        Sdata[5]    =   7.8e-21   ; 
edata[6]   =  5.2  ;        Sdata[6]     =   1.4e-20   ; 
edata[7]   =  8.    ;        Sdata[7]    =   2.5e-20   ; 
edata[8]   =  10.  ;        Sdata[8]    =   3.3e-20   ; 
edata[9]   =  13.  ;        Sdata[9]    =   4.5e-20   ; 
edata[10]  =  24.  ;        Sdata[10]   =   5.8e-20   ; 
edata[11]  =  64.  ;        Sdata[11]   =   6.2e-20   ; 
edata[12]  =  97.  ;        Sdata[12]   =   6.4e-20   ; 
edata[13]  =  147.  ;        Sdata[13]   =  6.6e-20    ; 
edata[14]  =  193.  ;        Sdata[14]   =  6.8e-20    ; 



 Emin = 0.005;             //eV
 Eth =  0.0;            //eV
 Emax = 40.005;         //eV
 Efc =  0.0;          //eV


//Ecoeff = (Omega_pe/dx)*(Omega_pe/dx)*(M_ions[4]+M_ions[13])/M_ions[13]/M_ions[4]/(0.5*T_e);

 Ecoeff = (dt/dr)*(dt/dr)*(mni_over_me+mn_over_me)/mn_over_me/mni_over_me/(0.5*T_e0);

 React_O2_n_detach.Emin = Emin*Ecoeff;
 React_O2_n_detach.Eth =  Eth*Ecoeff;
 React_O2_n_detach.Emax = Emax*Ecoeff;
// React_O2_n_detach.Wfc =  Efc*Ecoeff;            // square of e - O2+ relative velocity, corresponding to FC energy
 React_O2_n_detach.N = Npoints;
 React_O2_n_detach.Estep =  (React_O2_n_detach.Emax -  React_O2_n_detach.Emin)/Npoints;

  outputfile = fopen("xsct/O2_n_detach.dat", "w");

 for ( i = 0; i <= React_O2_n_detach.N; i++ )
 {
 E = i*React_O2_n_detach.Estep + React_O2_n_detach.Emin;

 tmp = 0.;

 if      (E/Ecoeff <= edata[0])  tmp = 0.0;
 else if (E/Ecoeff >= edata[14]) tmp = Sdata[14];
     else
   {
     int ind = 0;
    while (E/Ecoeff > edata[ind] ) ind++;
    assert( ind >= 1 && "calculated index cases out of bounds access to Sdata" );
    tmp = Sdata[ind - 1] + (Sdata[ind]- Sdata[ind-1])*(E/Ecoeff - edata[ind-1])/(edata[ind]- edata[ind-1]);
   }

 React_O2_n_detach.CS[i] = CScoeff*( tmp*1.e4 + 5.96*1.e-16/sqrt(E/Ecoeff)*0.17  )*coll_fac_ntrl_ntrl;   //coeff for low-energy part

    fprintf(outputfile, "%f  %15.7e \n", E/Ecoeff, React_O2_n_detach.CS[i]/CScoeff  );
   }
 fclose(outputfile); 
 
 
 /***********************************************************

        Mutual neutralization  // O- + O2+ = O + O2	
		Franz' analytic firmula
************************************************************/
 
 Emin = 0.002;             //eV
 Eth =  0.0;            //eV
 Emax = 16.002;         //eV
 Efc =  0.0;          //eV


//Ecoeff = (Omega_pe/dx)*(Omega_pe/dx)*(M_ions[4]+M_ions[6])/M_ions[6]/M_ions[4]/(0.5*T_e);

Ecoeff = (dt/dr)*(dt/dr)*(mni_over_me+mi_over_me)/mi_over_me/mni_over_me/(0.5*T_e0);
Ecoeff *= dt_ion*dt_ion;


 React_O2_Op_Om_ntrlzFr.Emin = Emin*Ecoeff;
 React_O2_Op_Om_ntrlzFr.Eth =  Eth*Ecoeff;
 React_O2_Op_Om_ntrlzFr.Emax = Emax*Ecoeff;
// React_O2_Op_Om_ntrlzFr.Wfc =  Efc*Ecoeff;            // square of e - O2+ relative velocity, corresponding to FC energy
 React_O2_Op_Om_ntrlzFr.N = Npoints;
 React_O2_Op_Om_ntrlzFr.Estep =  (React_O2_Op_Om_ntrlzFr.Emax -  React_O2_Op_Om_ntrlzFr.Emin)/Npoints;

  outputfile = fopen("xsct/O2_Op_Om_ntrlzFr.dat", "w");

 for ( i = 0; i <= React_O2_Op_Om_ntrlzFr.N; i++ )
 {
  E = i*React_O2_Op_Om_ntrlzFr.Estep + React_O2_Op_Om_ntrlzFr.Emin;

  React_O2_Op_Om_ntrlzFr.CS[i] = CScoeff*0.8*(1.0 + 2.85*Ecoeff/E)*1.e-14;  

    fprintf(outputfile, "%f  %15.7e \n", E/Ecoeff, React_O2_Op_Om_ntrlzFr.CS[i]/CScoeff  );
   }
 fclose(outputfile); 
}
