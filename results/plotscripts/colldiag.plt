nstep = time2d;

## seti is y, hubbel is x
seti = 3;
hubbel = 3;
xeich=1256
yeich=1024
x=hubbel*xeich;
y=seti*yeich;
print "size ".x.",".y
print "multiplot layout ".seti.",".hubbel

## input parameters for plot
set terminal pngcairo size x,y enhanced font 'PT Serif,30'
set termoption dash
set output '../../figs/'.twodrunid.'.collfigs.'.nstep.'.png'
set samples samplesiso
set isosamples samplesiso

## PLOT DEFINITIONS
load 'jet.pal'
set grid
unset key
set border linewidth 4.0
set tics nomirror out scale 1.5

## COLORS AND DASHES 
## define line styles using explicit rgbcolor names
set style line 1 lt 1 lc rgb "black" lw 2
set style line 2 lt 2 lc rgb "black" lw 2
set style line 3 lt 5 lc rgb "black" lw 2
set style line 4 lt 1 lc rgb "red" lw 2
set style line 5 lt 2 lc rgb "red" lw 2
set style line 6 lt 5 lc rgb "red" lw 2
set style line 7 lt 1 lc rgb "blue" lw 2
set style line 8 lt 2 lc rgb "blue" lw 2
set style line 9 lt 5 lc rgb "blue" lw 2
set style line 10 lt 1 lc rgb "cyan" lw 2
set style line 11 lt 2 lc rgb "cyan" lw 2
set style line 12 lt 5 lc rgb "cyan" lw 2

## label and axi
set multiplot layout seti,hubbel 
set xlabel 'z (cell)'
set ylabel 'r (cell)'
set logscale cb
set format cb "10^{%T}"

## PLOTS ITSELF
set title "coulomb coll at ".nstep
p '../out/coll_coulomb'.nstep.'.dat' matrix with image 
set title "electron neutral coll at ".nstep
p '../out/coll_el_ntrl'.nstep.'.dat' matrix with image 
set title "electron neutral excitation coll at ".nstep
p '../out/coll_el_ntrl_exc'.nstep.'.dat' matrix with image 
set title "collisional ionization at ".nstep
p '../out/coll_ioniz'.nstep.'.dat' matrix with image 
set title "ion neutral collisions at ".nstep
p '../out/coll_ion_ntrl'.nstep.'.dat' matrix with image 
set title "mean free path coll at ".nstep
p '../out/coll_mfp'.nstep.'.dat' matrix with image 
set title "nion neutral coll at ".nstep
p '../out/coll_nion_ntrl'.nstep.'.dat' matrix with image 
## RESET CB
unset logscale cb
if (debyediag == 1) {
	set title "local debye length from temp. (total energy) at ".nstep
	p '../out/e_debye_'.nstep.'.dat' matrix with image
}
set title "potential at ".nstep
p '../out/phi'.nstep.'.dat' matrix with image 

unset multiplot
