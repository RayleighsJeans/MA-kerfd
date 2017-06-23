function octave;
more off; warning off;
close all; clear all;

f_RF = 13.65e6;
source variables.tmp;

% PARTICLE NORM AND COLLFAC
particle_i = Ncell1/2;
particle_a = particle_i;
for i=1:nr
	particle_io = particle_i;
	particle_i = particle_io + Ncell1/2;
	particle_ao = particle_a;
	particle_a = particle_ao + particle_i;
end
particlenorm = particle_a/nr; 
dt_nion=dt_ion;

L_db02d = 7.43e2 * sqrt(Te02d / ne02d);
L_db01d = 7.43e2 * sqrt(Te01d / ne01d);
Omega_pe02d = 5.64e4 * sqrt(ne02d);
Omega_pe01d = 5.64e4 * sqrt(ne01d);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

collfac2d = pressure / (1.38e-23 * 300 * ne02d * 1e6) * (particlenorm / nn_pc);

%% physics stuff
elementary_chrg = 1.60219e-19;
% in Coulomb
atom_mass = 32.0;
% atom mass in au
mn_over_me = 1836.0 * atom_mass; 
% neutral mass over electron mass
mi_over_me = 1836.0 * atom_mass - 1.0;
% ion mass over electron mass
me_over_mi = 1.0 / mi_over_me;
% electron mass of ion mass
vac_permit = 8.854e-12;
% in A^2*s^4/kg*m^3, equals F/m
boltzman = 1.3806e-23;
% in J/K
v_te = dt / dr;

dr_0 = L_db02d*dr;
% cell width
dt_0 = dt/Omega_pe02d;
% time step length
dr1d_0 = L_db01d*dr;
% cell width 1D
explength1d = nz1d*dr1d_0;
explength = sizez2d*dr_0;
% experiment length in cm
expdiameter = nr*dr_0;
% experiment diameter in cm

n_RF2d = 1 / (f_RF * 0.2 / (Omega_pe02d));
% number of cycles for RF period 2D
n_RF1d = 1 / (f_RF * 0.2 / (Omega_pe01d));

%% AXIAL PROFILE APPROX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ndens1d = 1.25*collfac1d*ne01d;
ndens2d = 1.25*collfac2d*ne02d*fac1d2d;
suggestedfac = ndens2d/ndens1d;

%% REST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% electron velocity in dimensionless units
v_ti = v_te * sqrt(Ti_over_Te2d * me_over_mi) * dt_ion;
% ion velocity (thermal) dimensionless
v_tn = v_te * sqrt(Ti_over_Te2d * me_over_mi) * dt_ntrl;
% neutral thermal velocity dimensionless
% plasma frequency in Hz
% number of cycles for RF period 1D
Te02d_K = Te02d*elementary_chrg/boltzman;
% in K
T_i0 = Ti_over_Te2d*Te02d_K;
% ion temperature in K
ntr_fly = sqrt(3*boltzman*T_i0/(1.66e-27*32))*(0.2/Omega_pe02d)*(1/(0.5*L_db02d*1e-2))*dt_ion;
% neutral/ion flight width
debye_L1 = 1e2/sqrt((elementary_chrg*elementary_chrg)/(vac_permit*boltzman)*(ne02d*1e6)/(T_i0)*(1+Ti_over_Te2d));
% debye length test in cm
dr_0test1 = dr*debye_L1;
% cell width test in cm
debye_L2 = 1e2*sqrt(vac_permit*boltzman*Te02d_K/(ne02d*1e6*elementary_chrg*elementary_chrg));
% in cm
dr_0test2 = dr*debye_L2;
% cell width test in cm
Ndb = ne02d * 4/3*pi*L_db02d*L_db02d*L_db02d;
% number of particles in debye sphere
Ndb2 = ne02d * L_db01d;
% false number of particles in debye cell

%% SPEED OF SOUND FOR SPECIES ++ 1D DVT's %%%%%%%%%%%%%
%% SI SPEEDS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mi_over_me = 1832*atom_mass - 1;
mn_over_me = 1832*atom_mass;
v_te = (dt/dr)*dr_0/100./dt_0;
v_ti = v_te/dt_ion;
v_tni = v_te/dt_nion;
v_tn = v_te/dt_ntrl;
cs_e = (1/mi_over_me)**(1/2)*v_te;
cs_ntrl = dt_ntrl*cs_e;
cs_ion = dt_ion*cs_e;
cs_nion = dt_nion*cs_e;

%% FILES AND DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VELOCITIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1D
if (onedstuff == 1) 
  if (velzdiag == 1)
	  printf('>> radial velocities 1D\n');
	  for timecode=[timevec1d, time1d];
		  test=timecode+0;
		  if (test >= mintime1d) %mintime1d
			  if (test <= maxtime1d) %maxtime1d 
			  	printf('...%d',test);
		      file=strcat(onedfolder,'uex',num2str(timecode,'%08d'),'.dat'); uex=load(file);
		      file=strcat(onedfolder,'uey',num2str(timecode,'%08d'),'.dat'); uey=load(file);
		      file=strcat(onedfolder,'uO2px',num2str(timecode,'%08d'),'.dat'); uO2px=load(file);
		      file=strcat(onedfolder,'uO2py',num2str(timecode,'%08d'),'.dat'); uO2py=load(file);
		      file=strcat(onedfolder,'uOmx',num2str(timecode,'%08d'),'.dat'); uOmx=load(file);
		      file=strcat(onedfolder,'uOmy',num2str(timecode,'%08d'),'.dat'); uOmy=load(file);
	
		      uer=uO2pr=uOmr=zeros(length(uex(:,1)),2);
		      uer(:,1)=uex(:,1); uer(:,2)=0.5*sqrt(uey(:,2).^2+uex(:,2).^2);
		      file=strcat('transpose/uer',num2str(timecode,'%08d'),'.dat'); save("-text",file,"uer");
		      uO2pr(:,1)=uO2px(:,1); uO2pr(:,2)=0.5*sqrt(uO2py(:,2).^2+uO2px(:,2).^2);
		      file=strcat('transpose/uO2pr',num2str(timecode,'%08d'),'.dat'); save("-text",file,"uO2pr");
		      uOmr(:,1)=uOmx(:,1); uOmr(:,2)=0.5*sqrt(uOmy(:,2).^2+uOmx(:,2).^2);
		      file=strcat('transpose/uOmr',num2str(timecode,'%08d'),'.dat'); save("-text",file,"uOmr");
			  end
		  end
	  end
	  printf('\n');
  end
end
%% SAVING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear timevec1d timevec2d file uer uOmr uO2pr uex uOmx uO2px uOmy uey uO2py;
save -text 'data.dat' *
end 
