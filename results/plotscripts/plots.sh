#!/bin/bash 
#source variables.tmp

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
	echo ">> 1D/2D COMPARE"; gnuplot -e "${GnuVars}" compare1D2D.gplt;
	echo ">> PHYSFIGS"; gnuplot -e "${GnuVars}" physfigs.gplt;
	if [ "${colldiag}" = "1" ]; then echo ">> COLLDIAG"; gnuplot -e "${GnuVars}" colldiag.gplt; fi 
#	echo ">> XSECTIONS PLOTS"; gnuplot -e "${GnuVars}" xsections.gplt
	if [ "${velzdiag}" = "1" ]; then echo ">> VELZ COMPARE"; gnuplot -e "${GnuVars}" velzcompare.gplt; fi
fi
