#!/bin/bash

echo "PLOTS ################################################################"
echo "## PHYSCALC"
GnuVars+="nocalc='1'; "
gnuplot -e "${GnuVars}" calc.gplt

## plotting electron and positive ion density
## plotting density
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
