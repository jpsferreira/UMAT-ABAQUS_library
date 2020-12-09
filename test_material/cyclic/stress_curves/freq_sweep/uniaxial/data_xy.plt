set xlabel " strain"
set ylabel " stress"
m="./1.out"
set terminal postscript eps enhanced color font 'Helvetica,12'
set output 'laos.eps'
#set terminal png
#set output 'passive_vs_active_shear.png'
#set terminal postscript eps enhanced color font 'Helvetica,12'
#set output 'passive_vs_active_shear.eps'
#set nokey
set key bottom right
set grid
#set log y
#set log x
plot m using 2:3  with lines title 'laos'
#plot m using 1:3  with lines title 'passive',m using 1:4 with lines title 'active'


#set xtics nomirror
#set x2tics
#set autoscale xfix
#set autoscale x2fix
set terminal postscript eps enhanced color font 'Helvetica,12'
set output 'laos.eps'
set key bottom right
set grid
#
set multiplot
set size 1, 0.5
m="./1.out"

set origin 0.0,0.5
set xlabel "time"
set ylabel "stress"
set y2tics
set y2label "stretch"
#set log y
#set log x
plot m using 1:2 with lines title 'stretch' axis x1y2, \
     m using 1:3 with lines title 'stress' axis x1y1
set origin 0.0,0.0
set xlabel "stretch"
set ylabel "stress"
#set log y
#set log x
unset y2tics
unset y2label
plot m using 2:3 with lines title 'uniaxial'
unset multiplot
#plot m using 1:3  with lines title 'passive',m using 1:4 with lines title 'active'
#plot    m using 2:3  with lines title 'Lissajous'
