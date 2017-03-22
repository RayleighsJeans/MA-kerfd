#!/bin/bash
#source variables.tmp

##############################################################################
## TRANSPOSING ###############################################################
##############################################################################
## POTENTIAL AND DENSITY #####################################################
echo ">> transposing"

#echo ">> transposing potential and density: "
for time2d in ${timevec2d[@]}
do
	for name2d in ${names2d[@]}
	do
		 #echo -ne "... ${name2d}${time2d}.dat ";
		 head -${fullines} ../out/${name2d}${time2d}.dat | tail -${lines} > tmp.dat
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
		 }' tmp.dat  > transpose/trans-${name2d}${time2d}.dat
	 done
done

#echo ">> transposing velocities: "
for time2d in ${timevec2d[@]}
do
	for name2d in ${namev2d[@]}
	do
		for vsuf2d in ${velsuf2d[@]}
		do
			#echo -ne "... ${name2d}${vsuf2d}${time2d}.dat ";
			head -${fullines} ../out/${name2d}${vsuf2d}${time2d}.dat | tail -${lines} > tmp.dat
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
			}' tmp.dat  > transpose/trans-${name2d}${vsuf2d}${time2d}.dat
		done
	done
done

rm tmp.dat
echo ""
