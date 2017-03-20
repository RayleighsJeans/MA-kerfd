function octave_pac
more off; warning off;
close all; clear all;

%set(0,'defaultfigurevisible','off');
%set(0,'defaultaxesfontsize',24);
%set(0,'DefaultLinelinewidth',5);
%set(0,'Defaultaxeslinewidth',4);
%set(0,'defaulttextfontsize',24);
%set(0,'defaulttextfontweight','bold');
%axis('tick')

source variables.tmp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Omega_pe02d = 5.64e4 * sqrt(ne02d);
Omega_pe01d = 5.64e4 * sqrt(ne01d);
L_db02d = 7.43e2 * sqrt(Te02d / ne02d);
L_db01d = 7.43e2 * sqrt(Te01d / ne01d);
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
expdiameter = nr*dr_0;
% experiment diameter in cm
explength = sizez2d*dr_0;
% experiment length in cm
explength1d = nz1d*dr1d_0;

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
f_RF = 13.65e6;
% plasma frequency in Hz
n_RF2d = 1 / (f_RF * 0.2 / (Omega_pe02d));
% number of cycles for RF period 2D
n_RF1d = 1 / (f_RF * 0.2 / (Omega_pe01d));
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

%% FILES AND DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VELOCITIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1D
file=strcat(onedfolder,'uex',num2str(time1d,'%08d'),'.dat'); uex=load(file);
file=strcat(onedfolder,'uey',num2str(time1d,'%08d'),'.dat'); uey=load(file);
file=strcat(onedfolder,'uez',num2str(time1d,'%08d'),'.dat'); uez=load(file);

file=strcat(onedfolder,'uO2px',num2str(time1d,'%08d'),'.dat'); uO2px=load(file);
file=strcat(onedfolder,'uO2py',num2str(time1d,'%08d'),'.dat'); uO2py=load(file);
file=strcat(onedfolder,'uO2pz',num2str(time1d,'%08d'),'.dat'); uO2pz=load(file);

file=strcat(onedfolder,'uOmx',num2str(time1d,'%08d'),'.dat'); uOmx=load(file);
file=strcat(onedfolder,'uOmy',num2str(time1d,'%08d'),'.dat'); uOmy=load(file);
file=strcat(onedfolder,'uOmz',num2str(time1d,'%08d'),'.dat'); uOmz=load(file);

%% 2D
file='transpose2develr.dat'; evr=load(file); evelr=(1/lines)*sum(evr(:,:),2);
file='transpose2develz.dat'; evz=load(file); evelz=(1/lines)*sum(evz(:,:),2);
file='transpose2divelr.dat'; ivr=load(file); ivelr=(1/lines)*sum(ivr(:,:),2);
file='transpose2divelz.dat'; ivz=load(file); ivelz=(1/lines)*sum(ivz(:,:),2);
file='transpose2dnivelr.dat'; nivr=load(file); nivelr=(1/lines)*sum(nivr(:,:),2);
file='transpose2dnivelz.dat'; nivz=load(file); nivelz=(1/lines)*sum(nivz(:,:),2);
file='transpose2dnvelr.dat'; nvr=load(file); nvelr=(1/lines)*sum(nvr(:,:),2);
file='transpose2dnvelz.dat'; nvz=load(file); nvelz=(1/lines)*sum(nvz(:,:),2);

%% PROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uer=uO2pr=uOmr=zeros(length(uex(:,1)),2);
uer(:,1)=uex(:,1); uer(:,2)=sqrt(uey(:,2).^2+uex(:,2).^2);
uO2pr(:,1)=uO2px(:,1); uO2pr(:,2)=sqrt(uO2py(:,2).^2+uO2px(:,2).^2);
uOmr(:,1)=uOmx(:,1); uOmr(:,2)=sqrt(uOmy(:,2).^2+uOmx(:,2).^2);
length1d=uex(:,1).*L_db01d;
length2d=linspace(0.0,explength,length(evelr));

%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; xlabel('z/cm'); ylabel('v, dim.less');
		subplot(2,1,1);
		hold on;
		figtitle=strcat('electron velocities step-2D=',num2str(nstep2d),
								 ' nstep-1D=',num2str(nstep1d)); title(figtitle);
			plot(length1d,uez(:,2),'-.*k','linewidth',4);
			plot(length1d,uer(:,2),'--^k','linewidth',4);
			plot(length2d,evelr(:,1),'-.+r','linewidth',4);
			plot(length2d,evelz(:,1),'--pr','linewidth',4);
			legend('e velz 1D','e velr 1D','e velr 2D','e velz 2D'); hold off;
		subplot(2,1,2); hold on;
		figtitle=strcat('ion velocities step-2D=',num2str(nstep2d),
								 ' nstep-1D=',num2str(nstep1d)); title(figtitle);
			plot(length1d,uO2pz(:,2),'-.*k','linewidth',4);
			plot(length1d,uO2pr(:,2),'--^k','linewidth',4);
			plot(length2d,ivelr(:,1),'-.^+r','linewidth',4);
			plot(length2d,ivelz(:,1),'--pr','linewidth',4);
			legend('O2p velz 1D','O2p velr 1D','i velr 2D','i velz 2D'); hold off;
			hold off;
%		subplot(4,1,3); hold on;
%		figtitle=strcat('negative ion velocities step-2D=',num2str(nstep2d),
%								 ' nstep-1D=',num2str(nstep1d)); title(figtitle);
%			plot(length1d,uOmz(:,2),'-.*k','linewidth',4);
%			plot(length1d,uOmr(:,2),'--^k','linewidth',4);
%			plot(length2d,nivelr(:,1),'-.+r','linewidth',4);
%			plot(length2d,nivelz(:,1),'--pr','linewidth',4);
%			legend('Om velz 1D','Om velr 1D','ni velr 2D','ni velz 2D'); hold off;
%			hold off;
%		subplot(4,1,4); hold on;
%		figtitle=strcat('neutral velocities step-2D=',num2str(nstep2d)); 
%		title(figtitle);
%			plot(length2d,nvelr(:,1),'-.*r','linewidth',4);
%			plot(length2d,nvelz(:,1),'--pr','linewidth',4);
%			legend('n velr 2D','n velz 2D'); hold off;
%			hold off;
filename=strcat('figs/velz.1D-',num2str(ONEDRUNID),'.',num2str(nstep1d),
								'-2D.',num2str(TWODRUNID),'.',num2str(nstep2d),'.png');
print(filename,'-dpng','-S720,2160'); close;

keyboard;
end 
