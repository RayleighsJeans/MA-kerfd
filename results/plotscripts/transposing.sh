#!/bin/bash

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
