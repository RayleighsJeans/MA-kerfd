#!/bin/bash

timevec=( $(find | grep 'phi' | grep -Po '[0-9]+(?=.dat)' | sort -r | head -10 | sort -n) );
timesize=${#timevec[*]};
time=${timevec[$timesize-1]};
#time='04200000;'

# echo $time; echo $timesize;
# echo ${timevec[@]};
GnuVars+="time='${time}'; "
GnuVarsVec+="timevec= ' "; j="0";
	while [ $j -lt "$timesize" ]; do
		GnuVarsVec+="${timevec[$j]} ";
		j=$[$j+1];
	done
GnuVarsVec+="'; ";

GnuVars+="volt=400; "; #volts
GnuVars+="length=5.004; "; #cm
GnuVars+="Nx=426.0; "; #debye
GnuVars+="distnumber=320.0; " #standard
GnuVars+="ne0=5.0e9; "; #cm-3
GnuVars+="te0=5.0; "; #eV

join  Omin${time}.dat   Oms${time}.dat > 	 Om${time}.dat;
#join  Omin${time}.dat   Oms03700370.dat > 	 Om${time}.dat;
join   Tex${time}.dat   Tey${time}.dat >   Ter${time}.dat;
join   Ter${time}.dat   Tez${time}.dat >    Te${time}.dat;
join   Tix${time}.dat   Tiy${time}.dat >   Tir${time}.dat;
join   Tir${time}.dat   Tiz${time}.dat >    Ti${time}.dat;
join  Timx${time}.dat  Timy${time}.dat >  Timr${time}.dat;
join  Timr${time}.dat  Timz${time}.dat >   Tim${time}.dat;
join 	uez${time}.dat   uey${time}.dat > 	uer${time}.dat;
join uO2pz${time}.dat uO2py${time}.dat > uO2pr${time}.dat;
join  uOmz${time}.dat  uOmy${time}.dat >  uOmr${time}.dat;

gnuplot -e "${GnuVars} ${GnuVarsVec}" onedimensional.gplt;

exit
