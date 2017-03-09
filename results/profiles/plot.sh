#!/bin/bash
lines=60
ids=( "31306" "31307" "31309" "31314" "31315" );
prefix=( "n" "e" "i" );

GnuVars+="samplesiso=5000; "
GnuVars+="lines='${lines}'; "

j=0.0; k=0.0;

while [ j -le "${#prefix[@]}" ]; do
	GnuVars+="prefix[$j]='${prefix[$j]}'; "
	j=j+1;
done

while [ k -le "${#ids[@]}" ]; do
	GnuVars+="ids[$k]='${ids[$k]}'; "
	k=k+1;
done

echo ${GnuVars} > GnuVars.tmp

## TRANSPOSING ###############################################################
echo ">> >> loop for transpos. of files"
if ! [ "$1" -eq "$1" ] 2> /dev/null # checking if $1 is a number or not
	then 
		for i in "${ids[@]}"
		do
			for l in "${prefix[@]}"
			do
				cp -f ../../slurms/slurm-$i.out .
				cp -f ../w2d_${i}*/out/${l}_dens_00000015.dat $ldens$i.dat
				echo ">> doing transpositioning ${l}dens${i}.dat"
				## transpose
				head -${lines} ${l}dens${i}.dat > tmp.dat
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
				}' tmp.dat  > transposed/${l}dens${i}_trans.dat
			done
		done
fi

echo ">> >> doing plot/s"
gnuplot -e "${GnuVars}" gnuplot.gplt
