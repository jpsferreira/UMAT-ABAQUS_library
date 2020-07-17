#set xtics nomirror
#set x2tics
#set autoscale xfix
#set autoscale x2fix
set terminal postscript eps enhanced color font 'Helvetica,12'
set output 'monotonic.eps'
set key bottom right
set grid
#
set multiplot
set size 1, 0.5
m="./stress_curves/uniaxial.out"
p="./stress_curves/biaxial.out"
s="./stress_curves/shear.out"
ss="./stress_curves/s_shear.out"

set origin 0.0,0.5
set xlabel "Stretch"
set ylabel "Normal Stress"
plot m using 2:3 with lines title 'uniaxial', \
     p using 2:3 with lines title 'equibiaxial'
set origin 0.0,0.0
set xlabel "Amount of Shear"
set ylabel "Shear Stress"
plot s using 2:3 with lines title 'shear' , \
     ss using 2:3 with lines title 'simple shear' 
unset multiplot
#plot m using 1:3  with lines title 'passive',m using 1:4 with lines title 'active'
#plot    m using 2:3  with lines title 'Lissajous'
