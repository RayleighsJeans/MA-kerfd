#!/bin/bash 
# source variables.tmp
# source bash.vectors
# source Gnu.Vars

##############################################################################
## OCTAVE ####################################################################
echo ">> OCTAVE"; octave --silent --eval octave.m;

##############################################################################
## GIFS ######################################################################
if [ "1" = "${makegifs}" ]; then
	echo "GIFS ##################################################################" 
	gnuplot -e "${GnuVars} ${GnuVarsVec}" animation.gplt;
fi

##############################################################################
## FINALIZING AND PLOTS ######################################################
echo ""
if [ "1" = "${gnuplots}" ]; then
	echo ">> PLOTS"
	if [ ${onedstuff} = "1" ]; then echo ">> 1D/2D COMPARE"; gnuplot -e "${GnuVars} ${GnuVarsVec}" compare1D2D.gplt; fi
	echo ">> PHYSFIGS"; gnuplot -e "${GnuVars} ${GnuVarsVec}" physfigs.gplt;
	if [ "${colldiag}" = "1" ]; then echo ">> COLLDIAG"; gnuplot -e "${GnuVars} ${GnuVarsVec}" colldiag.gplt; fi 
  echo ">> XSECTIONS PLOTS"; gnuplot -e "${GnuVars} ${GnuVarsVec}" xsections.gplt
	if [ "${velzdiag}" = "1" ]; then
    if [ "${onedstuff}" = "1" ]; then
      echo ">> VELZ COMPARE"; gnuplot -e "${GnuVars} ${GnuVarsVec}" velzcompare.gplt;
    fi
  fi
fi
