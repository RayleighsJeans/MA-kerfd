#!/bin/bash

########################################################################
## TIMES AND LINES
########################################################################
copyonly=0; if [ "0" = "$copypnly" ]; then echo "## WILL ONLY COPY SHIT, THATS IT!"; fi
copyanything=1; if [ "0" = "$copyanything" ]; then echo "## WILL NOT COPY ANYTHING"; fi
copyh5s=1; if [ "0" = "$copyh5s" ]; then echo "## WILL NOT COPY *.H5's"; fi

outputnumber=25;
reduceoutput=1; if [ "1" = "$reduceoutput" ]; then echo "## WILL REDUCE OUTPUT TO ${outputnumber} lines"; fi
makegifs=1; if [ "0" = "$makegifs" ]; then echo "## WILL NOT MAKE GIFS"; fi
gnuplots=1; if [ "0" = "$gnuplots" ]; then echo "## WILL NOT MAKE GNUPLOTS"; fi
onedstuff=1; if [ "0" = "$onedstuff" ]; then echo "## WILL NOT DO ANY ONED STUFF"; fi

TWODRUNDIR=`exec pwd`
tmp=$(pwd | grep -Po '(_=?)+[0-9]+(?=.)')
TWODRUNID=${tmp##*[^0-9]}
TWODJOBNAME=$(grep -w "\#\# jobname\:" "../../slurms/slurm-${TWODRUNID}.out" | awk '{print($3)}')

if [ "$onedstuff" = "1" ]; then
  ## define which runs are to be compared
  ONEDRUNID=44387;
  ONEDRUNNAME=DoubleSpaceHP;
  ONEDRUNDIR=../../w1d_$ONEDRUNID.$ONEDRUNNAME;
fi

##############################################################################
## REMOVE SHIT AND STUFF #####################################################
rm -rf data.dat Gnu.Vars octave.vectors variables.tmp bash.vectors transpose/ 2>/dev/null;

##############################################################################
## COPY STUFF ################################################################
TWODBASEDIR=$TWODRUNDIR/../../../pic2d
TMPIC2DIR=$TWODBASEDIR/../tmpic2
TMPIC1DIR=$TWODBASEDIR/../tmpic1
SOLVER2DDIR=$TWODBASEDIR/../solver2d
if [ "$onedstuff" = "1" ]; then ONEDBASEDIR=$TWODBASEDIR/../pic1d; fi


if [ "1" = "$copyanything" ]; then
   echo ">> copy slurms";
   #head -3000 $SOLVER2DDIR/slurm-${TWODRUNID}.out > $TWODBASEDIR/../results/slurms/slurm-${TWODRUNID}.out 2>/dev/null;
   cp -fv $TWODBASEDIR/slurm-*.out $TWODBASEDIR/../results/slurms/ 2>/dev/null;
   if [ "$onedstuff" = "1" ]; then
     cp -fv $ONEDBASEDIR/slurm-*.out $TWODBASEDIR/../results/slurms/ 2>/dev/null;
   fi
   cp -fv $TMPIC2DIR/slurm-*.out $TWODBASEDIR/../results/slurms/ 2>/dev/null;
   cp -fv $TMPIC1DIR/slurm-*.out $TWODBASEDIR/../results/slurms/ 2>/dev/null;
   if [ "1" = "$copyh5s" ]; then
      echo ">> copy 2D particle *.h5's";
      cp -fv ../save/particle_backup0.h5 $TWODRUNDIR/../../particles/w2d_${TWODRUNID}-${TWODJOBNAME}_particles0.h5 2>/dev/null;
      cp -fv ../save/particle_backup1.h5 $TWODRUNDIR/../../particles/w2d_${TWODRUNID}-${TWODJOBNAME}_particles1.h5 2>/dev/null;
      cp -fv ../geom* $TWODRUNDIR/../../geometries/geom_${TWODRUNID}-${TWODJOBNAME}.txt 2>/dev/null
      cp -fv ../input.dat $TWODRUNDIR/../../inputs/2D-${TWODRUNID}-${TWODJOBNAME}.dat 2>/dev/null
      if [ "$onedstuff" = "1" ]; then
        cp -fv ${ONEDRUNDIR}/input.txt ${TWODRUNDIR}/../../inputs/1D-${ONEDRUNID}-${ONEDRUNNAME}.txt 2>/dev/null;
        cat ${ONEDRUNDIR}/subst.txt >> ${TWODRUNDIR}/../../inputs/1D-${ONEDRUNID}-${ONEDRUNNAME}.txt 2>/dev/null;
        echo ">> copy 1D particle *.dat's";
        cp -fv ${ONEDRUNDIR}/save/p0_0b.dat ${TWODRUNDIR}/../../particles/w1d_${ONEDRUNID}-${ONEDRUNNAME}_p0_0b.dat 2>/dev/null;
        cp -fv ${ONEDRUNDIR}/save/p0_1b.dat ${TWODRUNDIR}/../../particles/w1d_${ONEDRUNID}-${ONEDRUNNAME}_p0_1b.dat 2>/dev/null;
      fi
   fi
fi

########################################################################
## EVERYTHING ELSE IS DONE AFTER HERE ##################################
if [ "0" = "$copyonly" ]; then
  
  ########################################################################
  ## MISCELLANEOUS #######################################################
  axisoffset=10; fullines=60; lines=50;
  scale_arrow=10;
  names2d=( 'phi' 'e_dens_' 'i_dens_' 'ni_dens_' 'n_dens_');
  namev2d=( 'e_vel' 'i_vel' 'ni_vel' 'n_vel' ); velsuf2d=( 'r_' 'z_' 't_' );
  
  ##############################################################################
  ## PREPARATIONS ##############################################################
  mkdir -p figs 2>/dev/null;
  mkdir -p transpose 2>/dev/null;
  
  rm -rf transpose/* 2>/dev/null; # clear
  
  ##############################################################################
  ## shell input parameters ####################################################
  ## get 2D LAST NUMBER ########################################################
  if [ "$reduceoutput" = "1" ]; then
    echo ">> reduced timevec 2D";
    timevec2d=( $(ls ../out/phi* | grep -Po '[0-9]+(?=.dat)' | sort -r | head -${outputnumber} | sort -n) );
    echo ">> timevec2d: ${timevec2d[@]}";
  else
    timevec2d=( $(ls ../out/phi* | grep -Po '[0-9]+(?=.dat)') ); timesize2d=${#timevec2d[*]};
    echo ">> timevec2d: ${timevec2d[@]}";
  fi
  timesize2d=${#timevec2d[*]};
     echo ">> get 2D step number"
     if ! [ "$1" -eq "$1" ] 2> /dev/null; then 
        LAST2D=${timevec2d[$timesize2d-1]}; FIRST2D=${timevec2d[0]}; time2d=$LAST2D;
        nstep2d=`expr $time2d + 0`;  start2d=`expr $FIRST2D + 0`;
     else printf -v time2d "%08d" $nstep2d # given timestep
     fi
  
  ##############################################################################
  ## get 1D NUMBER closest to 2D number ########################################
  if [ "$onedstuff" = "1" ]; then
    cd $ONEDRUNDIR
    if [ "$reduceoutput" = "1" ]; then
      echo ">> reduced timevec 1D";
      timevec1d=( $(ls out/phi* | grep -Po '[0-9]+(?=.dat)' | sort -r | head -${outputnumber} | sort -n) );
      echo ">> timevec1d: ${timevec1d[@]}";
    else
      timevec1d=( $(ls out/ne* | grep -Po '[0-9]+(?=.dat)') );
      echo ">> timevec1d: ${timevec1d[@]}";
    fi
    timesize1d=${#timevec1d[*]};
    echo ">> get 1D step number"
    if ! [ "$1" -eq "$1" ] 2> /dev/null; then # checking if $1 is a number or not
       LAST1D=${timevec1d[$timesize1d-1]}; FIRST1D=${timevec1d[0]}; time1d=$LAST1D;
       nstep1d=`expr $time1d + 0`; start1d=`expr $FIRST1D + 0`;
        # echo ">> search 1D step number closest to 2D number"
        # if [ "${nstep2d}" != "${nstep1d}" ]; then
        #   i=0; base=100; newdiff=1; nstep1d=`expr ${timevec1d[$i]} + 0`;
        #   oldiff=`bc <<< "sqrt((${nstep2d} - ${nstep1d})^2)"`;
        #   echo -ne ">> >> iterating... "
        #   while [ "$oldiff" -ge "$newdiff" ]; do							
        #     i=$[$i+1];
        #     if [ $i -gt "1" ]; then oldiff=$newdiff; fi
        #     time1d=${timevec1d[$i]}; nstep1d=`expr ${time1d} + 0`;
        #     newdiff=`bc <<< "sqrt((${nstep2d} - ${nstep1d})^2)"`
        #     if [ $i -gt $base ]; then echo -ne "$newdiff ... "; base=$[$base+100]; fi
        #   done
        # fi
    else printf -v time1d "%08d" $1 # given timestep
    fi
  fi
  
  ##############################################################################
  ## IMPORTANTÈ ################################################################
  cd $TWODRUNDIR
  
  #time1d=00000002; nstep1d=2;
  #time2d=00000002; nstep2d=2;
  mintime1d=00000000; maxtime1d=1000000000000; #$nstep2d; #start2d;
  mintime2d=00000000; maxtime2d=1000000000000;
  
  ##############################################################################
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
  pressure=$(grep -w "pressure" "../input.dat" | awk '{print($2)}' | head -1) # pressure in pa
  sizez2d=$(grep -w "nz" "../input.dat" | awk '{print($2)}') # size in x direction
  ne02d=$(grep -w "n_e0" "../input.dat" | awk '{print($2)}') # ne standard electron density per cm³
  voltage2d=$(grep -w "Ua" "../input.dat" | awk '{print($2)}') # voltage2d at powered electrode
  Te02d=$(grep -w "T_e0" "../input.dat" | awk '{print($2)}') # electron temp in eV
  collfac2d=$(grep -oP '(?<=coll_fac_ntrl_ntrl=\ ).*' ../../slurms/slurm-$TWODRUNID.out | sort -n | tail -1)
  fac1d2d=$(grep -w "fac1d2d" "../input.dat" | awk '{print($2)}') # fac for 1D/2D neutral profile near axis
  Ti_over_Te2d=$(grep -w "Ti_over_Te" "../input.dat" | awk '{print($2)}') # electron temp in eV
  L_db02d=$(grep -oP '(?<=Debye-length=\ ).*' "../../slurms/slurm-$TWODRUNID.out" | sort -n | head -1) # debye length in cm
  taskpernode=$(exec cat ../../slurms/slurm-$TWODRUNID.out | grep "## SLURM TASKS PER NODE =" | grep -Po '[0-9]' | sort -n | tail -1)
  sizevec2d=( $(grep -oP "(?<=Global Parameters are\ ).*" "../../slurms/slurm-$TWODRUNID.out" | sort -n | head -1 ) );
  if [ -z $sizevec2d ]; then echo ">> sizevec is empty, use old paramters!";
    else
    echo ">> take global dimensions in cells from sizevec";
    sizer=${sizevec2d[0]}; sizez2d=${sizevec2d[1]};
    echo ">> global dimension nr=${sizer} and nz=${sizez2d}";
  fi
  samples=10000
  
  ##############################################################################
  ## 1D ########################################################################
  if [ "$onedstuff" = "1" ]; then 
  echo ">> get 1D parameters"
    nz1d=$(grep -w " nx" "../../slurms/slurm-$ONEDRUNID.out" | awk '{print($2)}' | head -1) # 1D debye length
    voltage1d=$(grep -w "#define SUBST_U0" "../../slurms/slurm-$ONEDRUNID.out" | grep -Po '[0-9]+(?=.*1.6/100.)' | tail -1) # applied voltag 1D 
    ne01d=$(grep -oP '(?<=#define SUBST_ne ).*' ../../slurms/slurm-$ONEDRUNID.out ) # standard density 1d
    Te01d=$(grep -oP '(?<=#define SUBST_Te ).*' ../../slurms/slurm-$ONEDRUNID.out ) # standard density 1d
    collfac1d=$(grep -oP '(?<=#define SUBST_ampl ).*' ../../slurms/slurm-$ONEDRUNID.out) # coll fac ntrls 1D from subst
    L_db01d=$(grep -w "Debye-length=" "../../slurms/slurm-$ONEDRUNID.out" | tail -1 | cut -d ' ' -f2-) ## debye length 1d cm
  fi
  
  ##############################################################################
  ## 2D FLAGS ##################################################################
  echo ">> get 2D flags"
  tempzdiag=$(grep -oP "(?<=DIAG_TEMPZ=).*" "../../slurms/slurm-$TWODRUNID.out" |  sort -n | tail -1) # tempz flag
  colldiag=$(grep -oP "(?<=DIAG_COLL=).*" "../../slurms/slurm-$TWODRUNID.out" | sort -n | tail -1) # colldiag flag
  velzdiag=$(grep -oP "(?<=DIAG_VELZ=).*" "../../slurms/slurm-$TWODRUNID.out" | sort -n | tail -1) # velz flag
  debyediag=$(grep -oP "(?<=DIAG_DEBYE_LENGTH=).*" "../../slurms/slurm-$TWODRUNID.out" | sort -n | tail -1) # debye flag
  ntrlzconst=$(grep -oP "(?<=NTRL_CONST=).*" "../../slurms/slurm-$TWODRUNID.out" | sort -n | tail -1) # ntrlz const
  negion=$(grep -oP "(?<=USE_FEEDGAS_OXYGEN=).*" "../../slurms/slurm-$TWODRUNID.out" | sort -n | tail -1) # neg ions 
  femsolver=$(grep -oP "(?<=USE_FEM_SOLVER=).*" "../../slurms/slurm-$TWODRUNID.out" | sort -n | tail -1) # neg ions 
  
  ##############################################################################
  ## 2D FLAG CHECK #############################################################
  true=$(echo "1")
  if [ "${velzdiag}" != "${true}" ]; then velzdiag=$(echo "0"); fi; ## VELZ
  if [ "${colldiag}" != "${true}" ]; then colldiag=$(echo "0"); fi ## COLLZ
  if [ "${negion}" != "${true}" ]; then negion=0; fi ## NEGATIVE IONS
  if [ "${tempzdiag}" != "${true}" ]; then tempzdiag=0; fi ## TEMPZ
  if [ "${debyediag}" != "${true}" ]; then debyediag=0; fi ## DEBYE
  if [ "${ntrlzconst}" != "${true}" ]; then ntrlzconst=1; fi ## NTRLZ STAY CONSTANT AND NO PLOT
  if [ "${femsolver}" != "${true}" ]; then femsolver=0; fi ## NTRLZ
  cd $TWODRUNDIR
  
  ##############################################################################
  ## GnuVars ###################################################################
  # 2D #########################################################################
  echo ">> 2D GnuVars"
  GnuVars+="dr='${dr2d}'; ";
  GnuVars+="dt='${dt2d}'; ";
  GnuVars+="atom_mass='${atom_mass2d}'; ";
  GnuVars+="collfac2d='${collfac2d}'; ";
  GnuVars+="fac1d2d='${fac1d2d}'; ";
  GnuVars+="L_db02d='${L_db02d}'; ";
  GnuVars+="Ti_over_Te2d='${Ti_over_Te2d}'; ";
  GnuVars+="Te02d='${Te02d}'; ";
  GnuVars+="voltage2d='${voltage2d}'; "
  GnuVars+="ne02d='${ne02d}'; "
  GnuVars+="sizez2d='${sizez2d}'; "
  GnuVars+="pressure='${pressure}'; "
  GnuVars+="dt_ntrl='${dt_ntrl}'; "
  GnuVars+="dt_ion='${dt_ion}'; "
  GnuVars+="nn_pc='${nn_pc}'; "
  GnuVars+="Ncell1='${Ncell1}'; "
  GnuVars+="navdt='${navdt}'; "
  GnuVars+="nr='${sizer}'; "
  GnuVars+="ntrlzconst='${ntrlzconst}'; "
  GnuVars+="debyediag='${debyediag}'; "
  GnuVars+="tempzdiag='${tempzdiag}'; "
  GnuVars+="negion='${negion}'; "
  GnuVars+="colldiag='${colldiag}'; "
  GnuVars+="velzdiag='${velzdiag}'; "
  GnuVars+="femsolver='${femsolver}'; "
  GnuVars+="onedstuff='${onedstuff}'; "
  GnuVars+="time2d='${time2d}'; " 
  GnuVars+="nstep2d='${nstep2d}'; "
  GnuVars+="twodfolder='${TWODRUNDIR}/../out/'; "
  GnuVars+="twodrunid='${TWODRUNID}'; "
  GnuVars+="scale_arrow='${scale_arrow}'; "
  GnuVars+="first2d='${FIRST2D}'; "
  GnuVars+="start2d='${start2d}'; "
  GnuVars+="last2d='${LAST2D}'; "
  GnuVars+="mintime2d='${mintime2d}'; "
  GnuVars+="maxtime2d='${maxtime2d}'; "
  
  ##############################################################################
  ## 1D ########################################################################
  if [ "$onedstuff" = "1" ]; then
    echo ">> 1D GnuVars"
    GnuVars+="time1d='${time1d}'; "
    GnuVars+="nstep1d='${nstep1d}'; "
    GnuVars+="onedrunid='${ONEDRUNID}'; "
    GnuVars+="onedfolder='${ONEDRUNDIR}/out/'; "
    GnuVars+="collfac1d='${collfac1d}'; "
    GnuVars+="ne01d='${ne01d}'; "
    GnuVars+="Te01d='${Te01d}'; "
    GnuVars+="nz1d='${nz1d}'; "
    GnuVars+="L_db01d='${L_db01d}'; "
    GnuVars+="voltage1d='${voltage1d}'; "
    GnuVars+="first1d='${FIRST1D}'; "
    GnuVars+="last1d='${LAST1D}'; "
    GnuVars+="mintime1d='${mintime1d}'; "
    GnuVars+="maxtime1d='${maxtime1d}'; "
  fi
  
  ##############################################################################
  ## MISC ######################################################################
  echo ">> MISC GnuVars"
  GnuVars+="lines='${lines}'; "
  GnuVars+="samplesiso='${samples}'; "
  GnuVars+="taskpernode='${taskpernode}'; "
  GnuVars+="nocalc='1'; "
  GnuVars+="copyonly='${copyonly}'; "
  GnuVars+="copyanything='${copyanything}'; "
  GnuVars+="makegifs='${makegifs}'; "
  GnuVars+="copyh5s='${copyh5s}'; "
  GnuVars+="gnuplots='${gnuplots}'; "
  GnuVars+="onedstuff='${onedstuff}'; "

  ##############################################################################
  ## DUMPING GnuVars to file ###################################################
  echo ">> dumping GnuVars to file";
      echo "GnuVars=${GnuVars};" >> Gnu.Vars;
  
  ##############################################################################
  ## ECHOING ###################################################################
  echo > variables.tmp
  ## 2D ########################################################################
  echo ""
  echo "2D STUFF ##############################################################"
  echo "time2d=${time2d}"; echo "time2d=${time2d};" >> variables.tmp;  
  echo "nstep2d=${nstep2d}"; echo "nstep2d=${nstep2d};" >> variables.tmp;
  echo "twodfolder=${TWODRUNDIR}/../out/"; echo "twodfolder='${TWODRUNDIR}/../out/';" >> variables.tmp; 
  echo "TWODRUNID=${TWODRUNID}"; echo "TWODRUNID=${TWODRUNID};" >> variables.tmp; 
  echo "dr2d=${dr2d}"; echo "dr=${dr2d};" >> variables.tmp; 
  echo "dt2d=${dt2d}"; echo "dt=${dt2d};" >> variables.tmp; 
  echo "atom_mass=${atom_mass2d}"; echo "atom_mass=${atom_mass2d};" >> variables.tmp; 
  echo "collfac2d=${collfac2d}"; echo "collfac2d=${collfac2d};" >> variables.tmp; 
  echo "fac1d2d=${fac1d2d}"; echo "fac1d2d=${fac1d2d};" >> variables.tmp; 
  echo "Ti_over_Te2d=${Ti_over_Te2d}"; echo "Ti_over_Te2d=${Ti_over_Te2d};" >> variables.tmp; 
  echo "Te02d=${Te02d}"; echo "Te02d=${Te02d};" >> variables.tmp; 
  echo "L_db02d=${L_db02d}"; echo "L_db02d=${L_db02d};" >> variables.tmp; 
  echo "voltage2d=${voltage2d}"; echo "voltage2d=${voltage2d};" >> variables.tmp; 
  echo "ne02d=${ne02d}"; echo "ne02d=${ne02d};" >> variables.tmp; 
  echo "sizez2d=${sizez2d}"; echo "sizez2d=${sizez2d};" >> variables.tmp; 
  echo "pressure=${pressure}"; echo "pressure=${pressure};" >> variables.tmp; 
  echo "dt_ntrl=${dt_ntrl}"; echo "dt_ntrl=${dt_ntrl};" >> variables.tmp; 
  echo "dt_ion=${dt_ion}"; echo "dt_ion=${dt_ion};" >> variables.tmp; 
  echo "nn_pc=${nn_pc}"; echo "nn_pc=${nn_pc};" >> variables.tmp; 
  echo "Ncell1=${Ncell1}"; echo "Ncell1=${Ncell1};" >> variables.tmp; 
  echo "navdt=${navdt}"; echo "navdt=${navdt};" >> variables.tmp; 
  echo "nr=${sizer}"; echo "nr=${sizer};" >> variables.tmp; 
  echo "scale_arrow=${scale_arrow}"; echo "scale_arrow=${scale_arrow};" >> variables.tmp; 
  echo "first2d=${FIRST2D}"; echo "first1d=${FIRST2D};" >> variables.tmp;
  echo "last2d=${LAST2D}"; echo "last2d=${LAST2D};" >> variables.tmp;
  echo "start2d=${start2d}"; echo "start2d=${start2d};" >> variables.tmp;
  echo "mintime2d=${mintime2d}"; echo "mintime2d=${mintime2d};" >> variables.tmp;
  echo "maxtime2d=${maxtime2d}"; echo "maxtime2d=${maxtime2d};" >> variables.tmp;
  echo ""
  
  ##############################################################################
  ## 1D ########################################################################
  if [ "$onedstuff" = "1" ]; then
    echo "1D STUFF ##############################################################"
    echo "time1d=${time1d}"; echo "time1d=${time1d};" >> variables.tmp;
    echo "nstep1d=${nstep1d}"; echo "nstep1d=${nstep1d};" >> variables.tmp;
    echo "onedfolder=${ONEDRUNDIR}/out/"; echo "onedfolder='${ONEDRUNDIR}/out/';" >> variables.tmp; 
    echo "ONEDRUNID=${ONEDRUNID}"; echo "ONEDRUNID=${ONEDRUNID};" >> variables.tmp; 
    echo "collfac1d=${collfac1d}"; echo "collfac1d=${collfac1d};" >> variables.tmp; 
    echo "ne01d=${ne01d}"; echo "ne01d=${ne01d};" >> variables.tmp; 
    echo "Te01d=${Te01d}"; echo "Te01d=${Te01d};" >> variables.tmp; 
    echo "voltage1d=${voltage1d}"; echo "voltage1d=${voltage1d};" >> variables.tmp; 
    echo "nz1d=${nz1d}"; echo "nz1d=${nz1d};" >> variables.tmp; 
    echo "L_db01d=${L_db01d}"; echo "L_db01d=${L_db01d};" >> variables.tmp; 
    echo "first1d=${FIRST1D}"; echo "first1d=${FIRST1D};" >> variables.tmp;
    echo "last1d=${LAST1D}"; echo "last1d=${LAST1D};" >> variables.tmp;
    echo "start1d=${start1d}"; echo "start1d=${start1d};" >> variables.tmp;
    echo "mintime1d=${mintime1d}"; echo "mintime1d=${mintime1d};" >> variables.tmp;
    echo "maxtime1d=${maxtime1d}"; echo "maxtime1d=${maxtime1d};" >> variables.tmp;
    echo "";  
  fi
  
  ##############################################################################
  ## FLAGS #####################################################################
  echo "FLAG STUFF ############################################################"
  echo "ntrlzconst=${ntrlzconst}"; echo "ntrlzconst=${ntrlzconst};" >> variables.tmp; 
  echo "debyediag=${debyediag}"; echo "debyediag=${debyediag};" >> variables.tmp; 
  echo "tempzdiag=${tempzdiag}"; echo "tempzdiag=${tempzdiag};" >> variables.tmp; 
  echo "negion=${negion}"; echo "negion=${negion};" >> variables.tmp; 
  echo "colldiag=${colldiag}"; echo "colldiag=${colldiag};" >> variables.tmp; 
  echo "velzdiag=${velzdiag}"; echo "velzdiag=${velzdiag};" >> variables.tmp; 
  echo "femsolver=${femsolver}"; echo "femsolver=${femsolver};" >> variables.tmp; 
  echo "onedstuff=${onedstuff}"; echo "onedstuff=${onedstuff};" >> variables.tmp; 
  echo "copyonly=${copyonly}"; echo "copyonly=${copyonly};" >> variables.tmp; 
  echo "copyh5s=${copyh5s}"; echo "copyh5s=${copyh5s};" >> variables.tmp; 
  echo "copyanything=${copyanything}";
       echo "copyanything=${copyanything};" >> variables.tmp; 
  echo "makegifs=${makegifs}"; echo "makegifs=${makegifs};" >> variables.tmp; 
  echo "gnuplots=${gnuplots}"; echo "gnuplots=${gnuplots};" >> variables.tmp; 
  echo "";  
  
  ##############################################################################
  ## MISC ######################################################################
  echo "MISC STUFF ############################################################"
  echo "taskpernode=${taskpernode}"; echo "taskpernode=${taskpernode};" >> variables.tmp; 
  echo "isosamples=${samples}"; echo "isosamples=${samples};" >> variables.tmp; 
  echo "lines=${lines}"; echo "lines=${lines};" >> variables.tmp; 
  echo "axisoffset=${axisoffset}"; echo "axisoffset=${axisoffset};" >> variables.tmp; 
  echo "fullines=${fullines}"; echo "fullines=${fullines};" >> variables.tmp; 
  echo "names2d=[ ${names2d[@]} ]";
  echo "namev2d=[ ${namev2d[@]} ]";
  echo "velsuf2d=[ ${velsuf2d[@]} ]";
  
  echo "names2d=( 'phi' 'e_dens_' 'i_dens_' 'ni_dens_' 'n_dens_' );" >> bash.vectors;
  echo "namev2d=( 'e_vel' 'i_vel' 'ni_vel' 'n_vel' );" >> bash.vectors;
  echo "velsuf2d=( 'r_' 'z_' 't_' );" >> bash.vectors;
  
  echo "names2d=[ 'phi' 'e_dens_' 'i_dens_' 'ni_dens_' 'n_dens_' ];" >> octave.vectors;
  echo "namev2d=[ 'e_vel' 'i_vel' 'ni_vel' 'n_vel' ];" >> octave.vectors;
  echo "velsuf2d=[ 'r_' 'z_' 't_' ];" >> octave.vectors;
  
  ##############################################################################
  ## TIMEVECTORS ###############################################################
  ## OCTAVE VECTORS ############################################################
  j="0";
  GnuVarsVec+="timevec2d= ' "; echo -ne "timevec2d=[ " >> octave.vectors;
  while [ $j -lt "$timesize2d" ]; do
    GnuVarsVec+="${timevec2d[$j]} ";
    echo -ne "${timevec2d[$j]} " >> octave.vectors;
    j=$[$j+1];
  done;
  GnuVarsVec+="'; "; echo "];" >> octave.vectors;
  echo "timesize2d=${timesize2d};"; GnuVarsVec+="timesize2d='${timesize2d}'; ";
  echo "timesize2d=${timesize2d};" >> variables.tmp; 
  
  if [ "$onedstuff" = "1" ]; then
    j="0";
    GnuVarsVec+="timevec1d= ' "; echo -ne "timevec1d=[ " >> octave.vectors;
    while [ $j -lt "$timesize1d" ]; do
      GnuVarsVec+="${timevec1d[$j]} ";
      echo -ne "${timevec1d[$j]} " >> octave.vectors;
      j=$[$j+5];
    done;
    GnuVarsVec+="'; "; echo "];" >> octave.vectors;
    echo "timesize1d=${timesize1d};"; GnuVarsVec+="timesize1d='${timesize1d}'; ";
    echo "timesize1d=${timesize1d};" >> variables.tmp; 
  fi
  echo ">> dumping GnuVarsVec";
      echo "GnuVarsVec=${GnuVarsVec};" >> Gnu.Vars;
  
##############################################################################
  ## BASH VECTORS ##############################################################
  
  j="0";
  echo -ne "timevec2d=( " >> bash.vectors;
  while [ $j -lt "$timesize2d" ]; do
    echo -ne "${timevec2d[$j]} " >> bash.vectors;
    j=$[$j+1];
  done;
  echo ");" >> bash.vectors;
  
  if [ "$onedstuff" = "1" ]; then
    j="0";
    echo -ne "timevec1d=( " >> bash.vectors
    while [ $j -lt "$timesize1d" ]; do
      echo -ne "${timevec1d[$j]} " >> bash.vectors;
      j=$[$j+5];
    done;
    echo ");" >> bash.vectors;
  fi
  
  ##############################################################################
  ## SUBSCRIPTS ################################################################
  echo ">> subscripting"
  . ./transpose.sh 
  . ./plots.sh 
  echo ">> done subscripting"
  
  cp -f data.dat figs/. 2>/dev/null;
  cp -f variables.tmp figs/. 2>/dev/null;
  mkdir -p ../../images/${TWODRUNID} 2>/dev/null;
  cp -rf figs/* ../../images/${TWODRUNID} 2>/dev/null;

fi
exit
