#!/bin/bash

TWODRUNDIR=`exec pwd`
tmp=$(pwd | grep -Po '(_=?)+[0-9]+(?=.)')
TWODRUNID=${tmp##*[^0-9]}
TWODJOBNAME=$(grep -w "\#\# jobname\:" "../../slurms/slurm-${TWODRUNID}.out" | awk '{print($3)}') # jobname

## COPY STUFF
TWODBASEDIR=$TWODRUNID/../../../tmpic2d
ONEDBASEDIR=$TWODBASEDIR/../pic1d

## cp $TWODBASEDIR/slurm-*.out $TWODBASEDIR/../results/slurms/ 2>/dev/null
## cp $ONEDBASEDIR/slurm-*.out $TWODBASEDIR/../results/slurms/ 2>/dev/null
## cp $TWODRUNDIR/../save/particle_backup0.h5 $TWODRUNDIR/../../particles/w2d_${TWODRUNID}-${TWODJOBNAME}_particles0.h5 2>/dev/null
## cp $TWODRUNDIR/../save/particle_backup1.h5 $TWODRUNDIR/../../particles/w2d_${TWODRUNID}-${TWODJOBNAME}_particles1.h5 2>/dev/null
## cp $TWODRUNDIR/../input.dat $TWODRUNDIR/../../inputs/2D-${TWODRUNID}-${TWODJOBNAME}.dat 2>/dev/null

## define which runs are to be compared
ONEDRUNID=31521;
ONEDRUNNAME=D288-6e4-31496;
ONEDRUNDIR=../../w1d_$ONEDRUNID.$ONEDRUNNAME

## shell input parameters #############################################
## get 2D LAST NUMBER #################################################
echo ">> get 2D step number"
if ! [ "$1" -eq "$1" ] 2> /dev/null # checking if $1 is a number or not
	then 
		LAST2D=`exec ls ../out/e_dens_* | grep -Po '[0-9]+(?=.dat)' | sort -n | tail -1`
		time2d=${LAST2D##*[^0-9]} # last timestep in folder out
		nstep2d=`expr $time2d + 0`; # timestep in nstep
		realtime2d=`awk -v "nstep=$nstep2d" '$2==nstep {print($1)}' ../trace.dat`; # clock time from nstep
	else
		printf -v time2d "%08d" $nstep2d # given timestep
fi
echo "time2d=${time2d}"
echo "nstep2d=${nstep2d}"
echo "realtime2d=${realtime2d}"

## get 1D NUMBER closest to 2D number ########################################
cd $ONEDRUNDIR
## echo ">> get 1D step number"
## if ! [ "$1" -eq "$1" ] 2> /dev/null # checking if $1 is a number or not
## 	then 
## 		LAST1D=`exec ls out/ | sed 's/\([0-9]\+\).*/\1/g' | sort -n | tail -1`
## 		time1d=${LAST1D##*[^0-9]} # last timestep in folder out
## 		nstep1d=`expr $time1d + 0`
## 		#echo "${nstep1d} == ${nstep2d} : " && bc <<< "${nstep1d} == ${nstep2d}"
## 		if [ "${nstep2d}" != "${nstep1d}" ]
## 			then
## 				i=1
## 				firstdiff=`bc <<< "sqrt((${nstep2d} - ${nstep1d})^2)"`
## 				oldiff=$firstdiff
## 				newdiff=1
## 				#echo "${firstdiff} >= ${newdiff} : " && bc <<< "${firstdiff} > ${newdiff}"
## 				while [ "$oldiff" -ge "$newdiff" ]
## 				do
## 					if [ $i -gt 1 ]; then
## 						oldiff=$newdiff
## 					fi
## 					OUT1D=`ls out/phi* | grep -Po '[0-9]+(?=.dat)' | sort -n | head -$i | tail -1`
## 					time1d=${OUT1D##*^[0-9]}
## 					nstep1d=`expr ${time1d} + 0`
## 					newdiff=`bc <<< "sqrt((${nstep2d} - ${nstep1d})^2)"`
## 					i=$[$i+1]
## 				done
## 	fi
## 	else
## 		printf -v time1d "%08d" $1 # given timestep
## fi
## 
## echo "time1d=${time1d}"

## IMPORTANTÈ ################################################################
cd $TWODRUNDIR

##############################################################################
## parameters from input.dat
## 2D ########################################################################
echo ">> get 2D parameters"
dr2d=$(grep -w "dr" "../input.dat" | awk '{print($2)}') # space sampling
dt2d=$(grep -w "dt" "../input.dat" | awk '{print($2)}') # time sampling
atom_mass2d=$(grep -w "atom_mass" "../input.dat" | awk '{print($2)}') # atom mass 
sizer=$(grep -w "nr" "../input.dat" | awk '{print($2)}') # size in y direction
navdt=$(grep -w "nav_dt" "../input.dat" | awk '{print($2)}') # size in y direction
Ncell1=$(grep -w "Ncell1" "../input.dat" | awk '{print($2)}') # Ncell1
nn_pc=$(grep -w "nn_pc" "../input.dat" | awk '{print($2)}') # nn_pc
dt_ion=$(grep -w "dt_ion" "../input.dat" | awk '{print($2)}') # dt_ion
dt_ntrl=$(grep -w "dt_ntrl" "../input.dat" | awk '{print($2)}') # dt_ntrl
dt_nion=$(grep -w "dt_nion" "../input.dat" | awk '{print($2)}') # dt_nion
pressure=$(grep -w "pressure" "../input.dat" | awk '{print($2)}' | head -1) # pressure in pa
sizez2d=$(grep -w "nz" "../input.dat" | awk '{print($2)}') # size in x direction
ne02d=$(grep -w "n_e0" "../input.dat" | awk '{print($2)}') # ne standard electron density per cm³
voltage2d=$(grep -w "Ua" "../input.dat" | awk '{print($2)}') # voltage2d at powered electrode
Te02d=$(grep -w "T_e0" "../input.dat" | awk '{print($2)}') # electron temp in eV
collfac2d=$(grep -oP '(?<=coll_fac_ntrl_ntrl=).*' ../../slurms/slurm-$TWODRUNID.out | sort -n | tail -1)
fac1d2d=$(grep -w "fac1d2d" "../input.dat" | awk '{print($2)}') # fac for 1D/2D neutral profile near axis
particlenorm=$(grep -oP '(?<=>> Density scaling neutrals: For particle norm \(nn_pc\) is ).*' "../../slurms/slurm-$TWODRUNID.out") # 2D particle norm
Ti_over_Te2d=$(grep -w "Ti_over_Te" "../input.dat" | awk '{print($2)}') # electron temp in eV
L_db02d=$(grep -oP '(?<=L_db0=).*.(?<=e-02)' "../../slurms/slurm-$TWODRUNID.out") # debye length in cm
taskpernode=$(exec cat ../../slurms/slurm-$TWODRUNID.out | grep "## SLURM TASKS PER NODE =" | grep -Po '[0-9]' | sort -n | tail -1)
samples=10000
## 1D ########################################################################
echo ">> get 1D parameters"
nz1d=$(grep -w " nx" "../../slurms/slurm-$ONEDRUNID.out" | awk '{print($2)}' | head -1) # 1D debye length
voltage1d=$(grep -w "#define SUBST_U0" "../../slurms/slurm-$ONEDRUNID.out" | grep -Po '[0-9]+(?=.*1.6/100.)' | tail -1) # applied voltag 1D 
ne01d=$(grep -oP '(?<=#define SUBST_ne ).*' ../../slurms/slurm-$ONEDRUNID.out ) # standard density 1d
Te01d=$(grep -oP '(?<=#define SUBST_Te ).*' ../../slurms/slurm-$ONEDRUNID.out ) # standard density 1d
collfac1d=$(grep -oP '(?<=#define SUBST_ampl ).*' ../../slurms/slurm-$ONEDRUNID.out) # coll fac ntrls 1D from subst
L_db01d=$(grep -w "Debye-length=" "../../slurms/slurm-$ONEDRUNID.out" | tail -1 | cut -d ' ' -f2-) ## debye length 1d cm
## 2D FLAGS ##################################################################
echo ">> get 2D flags"
rankdebug=$(grep -w "DEBUG_PARTICLE_CELLS_AND_RANKS" "../../slurms/slurm-$TWODRUNID.out" | awk '{print($2)}'| sort -n | tail -1) # if rankdebug flag is set
tempzdiag=$(grep -w "DIAG_TEMPZ" "../../slurms/slurm-$TWODRUNID.out" | awk '{print($2)}'| sort -n | tail -1) # if tempz flag is set
colldiag=$(grep -w "DIAG_COLL" "../../slurms/slurm-$TWODRUNID.out" | awk '{print($2)}'| sort -n | tail -1) # if colldiag flag is set
velzdiag=$(grep -w "DIAG_VELZ" "../../slurms/slurm-$TWODRUNID.out" | awk '{print($2)}'| sort -n | tail -1) # if velz flag is set
debyediag=$(grep -w "DIAG_DEBYE_LENGTH" "../../slurms/slurm-$TWODRUNID.out" | awk '{print($2)}'| sort -n | tail -1) # if debye flag is set
ntrlzconst=$(grep -w "NTRL_CONST" "../../slurms/slurm-$TWODRUNID.out" | awk '{print($2)}'| sort -n | tail -1) # if ntrlz const
reactdebug=$(grep -w "DEBUG_REACT" "../../slurms/slurm-$TWODRUNID.out" | awk '{print($2)}'| sort -n | tail -1) # if reaction debug
negion=$(grep -w "USE_NEGATIVE_IONS" "../../slurms/slurm-$TWODRUNID.out" | awk '{print($2)}'| sort -n | tail -1) # if negative ions 
ranknmb=$(cat ../../slurms/slurm-$TWODRUNID.out | grep "Rank #" | grep -Po '[0-9]+(?= is here and ready!)' | sort -n | tail -1) # who many ranks are used

##############################################################################
## 2D FLAG CHECK #############################################################
true=$(echo "1")
if [ "${velzdiag}" != "${true}" ] ## VELZ
then
	velzdiag=$(echo "0")
fi
if [ "${colldiag}" != "${true}" ] ## COLLZ
then
	colldiag=$(echo "0")
fi
if [ "${rankdebug}" != "${true}" ] ## RANKS
then
	rankdebug=0
fi
if [ "${rankdebug}" -eq "${true}" ] ## RAKDEBUG & CELLS
then
	nstepcells=`ls ../ | grep "neutrals" | grep -Po '[0-9]+(?=.cells.rank)' | sort -n | tail -1`
fi
if [ "${negion}" != "${true}" ] ## NEGATIVE IONS
then
	negion=0
fi
if [ "${tempzdiag}" != "${true}" ] ## TEMPZ
then
	tempzdiag=0
fi
if [ "${debyediag}" != "${true}" ] ## DEBYE
then
	debydiag=0
fi
if [ "${ntrlzconst}" != "${true}" ] ## NTRLZ
then
	ntrlzconst=0
fi
if [ "${reactdebug}" != "${true}" ] ## NTRLZ
then
	reactdebug=0
fi

cd $TWODRUNDIR
########################################################################
## TIMES AND LINES
########################################################################
lines=60
## time1d=00075000
## nstep1d=75000
## time2d=00016516
## nstep2d=16516

##############################################################################
## GnuVars
# 2D #########################################################################
echo ">> 2D GnuVars"
GnuVars+="dr='${dr2d}'; "; export dr;
GnuVars+="dt='${dt2d}'; "; export dt;
GnuVars+="atom_mass='${atom_mass2d}'; "; export atom_mass;
GnuVars+="particlenorm='${particlenorm}'; "; export particlenorm;
GnuVars+="collfac2d='${collfac2d}'; "export collfac2d;
GnuVars+="fac1d2d='${fac1d2d}'; "; export fac1d2d;
GnuVars+="L_db02d='${L_db02d}'; "; export L_db02d;
GnuVars+="Ti_over_Te2d='${Ti_over_Te2d}'; "; export Ti_over_Te2d;
GnuVars+="Te02d='${Te02d}'; "; export Te02d;
GnuVars+="voltage2d='${voltage2d}'; "
GnuVars+="ne02d='${ne02d}'; "
GnuVars+="sizez2d='${sizez2d}'; "
GnuVars+="pressure='${pressure}'; "
GnuVars+="dt_ntrl='${dt_ntrl}'; "
GnuVars+="dt_ion='${dt_ion}'; "
GnuVars+="dt_nion='${dt_nion}'; "
GnuVars+="nn_pc='${nn_pc}'; "
GnuVars+="Ncell1='${Ncell1}'; "
GnuVars+="navdt='${navdt}'; "
GnuVars+="nr='${sizer}'; "
GnuVars+="ntrlzconst='${ntrlzconst}'; "
GnuVars+="debyediag='${debyediag}'; "
GnuVars+="rank='${ranknmb}'; "
GnuVars+="tempzdiag='${tempzdiag}'; "
GnuVars+="negion='${negion}'; "
GnuVars+="rankdebug='${rankdebug}'; "
GnuVars+="nstepcells='${nstepcells}'; "
GnuVars+="colldiag='${colldiag}'; "
GnuVars+="velzdiag='${velzdiag}'; "
GnuVars+="reactdebug='${reactdebug}'; "
GnuVars+="realtime2d='${realtime2d}'; "
GnuVars+="time2d='${time2d}'; " # gnuplot variable of time 
GnuVars+="nstep2d='${nstep2d}'; "
GnuVars+="twodfolder='${TWODRUNDIR}/../out/'; "
GnuVars+="twodrunid='${TWODRUNID}'; "
## 1D ########################################################################
echo ">> 1D GnuVars"
GnuVars+="timeX1d='${time1d}'; " # gnuplot variable of time
GnuVars+="time1d='${nstep1d}'; "
GnuVars+="onedrunid='${ONEDRUNID}'; "
GnuVars+="onedfolder='${ONEDRUNDIR}/out/'; "
GnuVars+="collfac1d='${collfac1d}'; "
GnuVars+="ne01d='${ne01d}'; "
GnuVars+="Te01d='${Te01d}'; "
GnuVars+="nz1d='${nz1d}'; "
GnuVars+="L_db01d='${L_db01d}'; "
GnuVars+="voltage1d='${voltage1d}'; "
## MISC ######################################################################
echo ">> MISC GnuVars"
GnuVars+="lines='${lines}'; "
GnuVars+="samplesiso='${samples}'; "
GnuVars+="taskpernode='${taskpernode}'; "
GnuVars+="nocalc='1'; "

##############################################################################
## ECHOING
##############################################################################
## 2D ########################################################################
echo ""
echo "2D STUFF ##############################################################"
echo "twodfolder=${TWODRUNDIR}/../out/"
echo "TWODRUNID=${TWODRUNID}"
echo "dr=${dr2d}"
echo "dt=${dt2d}"
echo "atom_mass=${atom_mass2d}"
echo "particlenorm=${particlenorm}"
echo "collfac2d=${collfac2d}"
echo "fac1d2d=${fac1d2d}"
echo "Ti_over_Te2d=${Ti_over_Te2d}"
echo "Te0=${Te02d}"
echo "L_db02d=${L_db02d}"
echo "voltage2d=${voltage2d}"
echo "ne02d=${ne02d}"
echo "sizez2d=${sizez2d}"
echo "pressure=${pressure}"
echo "dt_ntrl=${dt_ntrl}"
echo "dt_ion=${dt_ion}"
echo "dt_nion=${dt_nion}"
echo "nn_pc=${nn_pc}"
echo "Ncell1=${Ncell1}"
echo "navdt=${navdt}"
echo "nr=${sizer}"
echo "time2d=${time2d}"
echo ""
## 1D ########################################################################
echo "1D STUFF ##############################################################"
echo "onedfolder=${ONEDRUNDIR}/out/"
echo "ONEDRUNID=${ONEDRUNID}"
echo "collfac1d=${collfac1d}"
echo "ne01d=${ne01d}"
echo "Te01d=${Te01d}"
echo "voltage1d=${voltage1d}"
echo "nz1d=${nz1d}"
echo "L_db01d=${L_db01d}"
echo ""
## FLAGS #####################################################################
echo "FLAG STUFF ############################################################"
echo "ntrlzconst=${ntrlzconst}"
echo "debyediag=${debyediag}"
echo "tempzdiag=${tempzdiag}"
echo "negion=${negion}"
echo "rankdebug=${rankdebug}"
echo "nstepcells=${nstepcells}" # last step in cells debug
echo "colldiag=${colldiag}"
echo "velzdiag=${velzdiag}"
echo "reactdebug=${reactdebug}"
echo ""
## MISC ######################################################################
echo "MISC STUFF ############################################################"
echo "ranknmb=${ranknmb}"
echo "taskpernode=${taskpernode}"
echo "isosamples=${samples}"
echo "lines=${lines}"

echo "TRANSPOSING ###########################################################"
## TRANSPOSING ###############################################################
##############################################################################
## POTENTIAL AND DENSITY #####################################################
## transpose potential
head -${lines} $TWODRUNDIR/../out/phi${time2d}.dat > tmp.dat
awk '
{
	for (i=1; i<=NF; i++){
		a[NR,i]=$i
	}
}
NF>p { p = NF }
END {
		for(j=1; j<=p; j++) {
			str=a[1,j]
			for(i=2; i<=NR; i++){
				str=str" "a[i,j];
			}
			print str
		}
}' tmp.dat  > transpose2dpot.dat
## transpose ion density
head -${lines} $TWODRUNDIR/../out/i_dens_${time2d}.dat > tmp.dat
awk '
{
	for (i=1; i<=NF; i++){
		a[NR,i]=$i
	}
}
NF>p { p = NF }
END {
		for(j=1; j<=p; j++) {
			str=a[1,j]
			for(i=2; i<=NR; i++){
				str=str" "a[i,j];
			}
			print str
		}
}' tmp.dat > transpose2didens.dat
## transpose neural density
head -${lines} $TWODRUNDIR/../out/n_dens_${time2d}.dat > tmp.dat
awk '
{
	for (i=1; i<=NF; i++){
		a[NR,i]=$i
	}
}
NF>p { p = NF }
END {
		for(j=1; j<=p; j++) {
			str=a[1,j]
			for(i=2; i<=NR; i++){
				str=str" "a[i,j];
			}
			print str
		}
}' tmp.dat > transpose2dnndens.dat
## transpose nion density
head -${lines} $TWODRUNDIR/../out/ni_dens_${time2d}.dat > tmp.dat
awk '
{
	for (i=1; i<=NF; i++){
		a[NR,i]=$i
	}
}
NF>p { p = NF }
END {
		for(j=1; j<=p; j++) {
			str=a[1,j]
			for(i=2; i<=NR; i++){
				str=str" "a[i,j];
			}
			print str
		}
}' tmp.dat > transpose2dnidens.dat
## transpose electron density
head -${lines} $TWODRUNDIR/../out/e_dens_${time2d}.dat > tmp.dat
awk '
{
	for (i=1; i<=NF; i++){
		a[NR,i]=$i
	}
}
NF>p { p = NF }
END {
		for(j=1; j<=p; j++) {
			str=a[1,j]
			for(i=2; i<=NR; i++){
				str=str" "a[i,j];
			}
			print str
		}
}' tmp.dat > transpose2dedens.dat

##############################################################################
## VELOCITIES ################################################################
## transpose electron r velocity
head -${lines} $TWODRUNDIR/../out/e_velr_${time2d}.dat > tmp.dat
awk '
{
	for (i=1; i<=NF; i++){
		a[NR,i]=$i
	}
}
NF>p { p = NF }
END {
		for(j=1; j<=p; j++) {
			str=a[1,j]
			for(i=2; i<=NR; i++){
				str=str" "a[i,j];
			}
			print str
		}
}' tmp.dat > transpose2develr.dat
## transpose electron z velocity
head -${lines} $TWODRUNDIR/../out/e_velz_${time2d}.dat > tmp.dat
awk '
{
	for (i=1; i<=NF; i++){
		a[NR,i]=$i
	}
}
NF>p { p = NF }
END {
		for(j=1; j<=p; j++) {
			str=a[1,j]
			for(i=2; i<=NR; i++){
				str=str" "a[i,j];
			}
			print str
		}
}' tmp.dat > transpose2develz.dat
## transpose electron thermal velocity
head -${lines} $TWODRUNDIR/../out/e_velt_${time2d}.dat > tmp.dat
awk '
{
	for (i=1; i<=NF; i++){
		a[NR,i]=$i
	}
}
NF>p { p = NF }
END {
		for(j=1; j<=p; j++) {
			str=a[1,j]
			for(i=2; i<=NR; i++){
				str=str" "a[i,j];
			}
			print str
		}
}' tmp.dat > transpose2develt.dat
## transpose ion r velocity
head -${lines} $TWODRUNDIR/../out/i_velr_${time2d}.dat > tmp.dat
awk '
{
	for (i=1; i<=NF; i++){
		a[NR,i]=$i
	}
}
NF>p { p = NF }
END {
		for(j=1; j<=p; j++) {
			str=a[1,j]
			for(i=2; i<=NR; i++){
				str=str" "a[i,j];
			}
			print str
		}
}' tmp.dat > transpose2divelr.dat
## transpose ion z velocity
head -${lines} $TWODRUNDIR/../out/i_velz_${time2d}.dat > tmp.dat
awk '
{
	for (i=1; i<=NF; i++){
		a[NR,i]=$i
	}
}
NF>p { p = NF }
END {
		for(j=1; j<=p; j++) {
			str=a[1,j]
			for(i=2; i<=NR; i++){
				str=str" "a[i,j];
			}
			print str
		}
}' tmp.dat > transpose2divelz.dat
## transpose ion themal velocity
head -${lines} $TWODRUNDIR/../out/i_velt_${time2d}.dat > tmp.dat
awk '
{
	for (i=1; i<=NF; i++){
		a[NR,i]=$i
	}
}
NF>p { p = NF }
END {
		for(j=1; j<=p; j++) {
			str=a[1,j]
			for(i=2; i<=NR; i++){
				str=str" "a[i,j];
			}
			print str
		}
}' tmp.dat > transpose2divelt.dat
## transpose nion r velocity
head -${lines} $TWODRUNDIR/../out/ni_velr_${time2d}.dat > tmp.dat
awk '
{
	for (i=1; i<=NF; i++){
		a[NR,i]=$i
	}
}
NF>p { p = NF }
END {
		for(j=1; j<=p; j++) {
			str=a[1,j]
			for(i=2; i<=NR; i++){
				str=str" "a[i,j];
			}
			print str
		}
}' tmp.dat > transpose2dnivelr.dat
## transpose electron z velocity
head -${lines} $TWODRUNDIR/../out/ni_velz_${time2d}.dat > tmp.dat
awk '
{
	for (i=1; i<=NF; i++){
		a[NR,i]=$i
	}
}
NF>p { p = NF }
END {
		for(j=1; j<=p; j++) {
			str=a[1,j]
			for(i=2; i<=NR; i++){
				str=str" "a[i,j];
			}
			print str
		}
}' tmp.dat > transpose2dnivelz.dat
## transpose nion thermal velocity
head -${lines} $TWODRUNDIR/../out/ni_velt_${time2d}.dat > tmp.dat
awk '
{
	for (i=1; i<=NF; i++){
		a[NR,i]=$i
	}
}
NF>p { p = NF }
END {
		for(j=1; j<=p; j++) {
			str=a[1,j]
			for(i=2; i<=NR; i++){
				str=str" "a[i,j];
			}
			print str
		}
}' tmp.dat > transpose2dnivelt.dat
## transpose neutral r velocity
head -${lines} $TWODRUNDIR/../out/n_velr_${time2d}.dat > tmp.dat
awk '
{
	for (i=1; i<=NF; i++){
		a[NR,i]=$i
	}
}
NF>p { p = NF }
END {
		for(j=1; j<=p; j++) {
			str=a[1,j]
			for(i=2; i<=NR; i++){
				str=str" "a[i,j];
			}
			print str
		}
}' tmp.dat > transpose2dnvelr.dat
## transpose neutral z velocity
head -${lines} $TWODRUNDIR/../out/n_velz_${time2d}.dat > tmp.dat
awk '
{
	for (i=1; i<=NF; i++){
		a[NR,i]=$i
	}
}
NF>p { p = NF }
END {
		for(j=1; j<=p; j++) {
			str=a[1,j]
			for(i=2; i<=NR; i++){
				str=str" "a[i,j];
			}
			print str
		}
}' tmp.dat > transpose2dnvelz.dat
## transpose neutral thermal velocity
head -${lines} $TWODRUNDIR/../out/n_velt_${time2d}.dat > tmp.dat
awk '
{
	for (i=1; i<=NF; i++){
		a[NR,i]=$i
	}
}
NF>p { p = NF }
END {
		for(j=1; j<=p; j++) {
			str=a[1,j]
			for(i=2; i<=NR; i++){
				str=str" "a[i,j];
			}
			print str
		}
}' tmp.dat > transpose2dnvelt.dat


## ANIMATION #################################################################
echo "GIFS ##################################################################" 
TIMEVEC=`exec ls ../out/e_dens_* | grep -Po '[0-9]+(?=.dat)' | sort -n`

j="0"

GnuVars+="timevec2d= ' ";
while [ $j -lt "${#TIMEVEC[@]}" ]; do
	GnuVars+="${TIMEVEC[$j]} ";
	j=$[$j+1];
done
GnuVars+="'; "
framenumber="${#TIMEVEC[@]}"
GnuVars+="framenumber='${framenumber}'; "

## making gifs
echo "## ANIMATION EDENS/IDENS/NIDENS"
gnuplot -e "${GnuVars}" animation.gplt

##############################################################################
## FINALIZING AND PLOTS ######################################################
echo ""
echo "PLOTS #################################################################"
##echo "## PHYSCALC"
##gnuplot -e "${GnuVars}" calc.gplt
## plotting electron and positive ion density
echo "## 1D/2D COMPARE"
gnuplot -e "${GnuVars}" compare1D2D.gplt
if [ "${rankdebug}" = "1" ]
	then
		echo "## CELL AND RANK DEBUG"
		gnuplot -e "${GnuVars}" cell_and_ranks.gplt
fi
echo "## PHYSFIGS"
gnuplot -e "${GnuVars}" physfigs_nstep.gplt
if [ "${colldiag}" = "1" ]
	then
		echo "## COLLDIAG"
		gnuplot -e "${GnuVars}" colldiag.plt
fi
gnuplot -e "${GnuVars}" xsections.gplt
if [ "${velzdiag}" = "1" ]
then
	echo "## VELZ COMPARE"
	gnuplot -e "${GnuVars}" velzcompare.gplt
fi

## DUMP ######################################################################
echo ${GnuVars} > GnuVars.tmp

rm GnuVars.tmp 2>/dev/null
rm tmp.dat 2>/dev/null
rm transpose* 2>/dev/null

exit
