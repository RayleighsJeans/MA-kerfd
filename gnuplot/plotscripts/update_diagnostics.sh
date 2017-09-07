#!/bin/bash

BASEDIR='/home/ph123693/Documents/2D/results/'
TWODRUNIDVEC=('44254' '44256' '44258' '44326' '44329' '44330' '44327' '44332' '44333');
cd ${BASEDIR}; 

for twodrunid in ${TWODRUNIDVEC[@]}
do
  echo ">> updating diag of ${twodrunid}"; 
  cd w2d_${twodrunid}*/gnuplot/;
  bash script_all.sh;
  cd ${BASEDIR};
done
