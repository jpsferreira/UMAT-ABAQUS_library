#set xtics nomirror
#set x2tics
#set autoscale xfix
#set autoscale x2fix
set terminal postscript eps enhanced color font 'Helvetica,12'
set output 'cyclic_freq.eps'
set key bottom right
set grid
#
set multiplot
set size 1, 0.5
m="./stress_curves/freq_sweep/uniaxial.out"
p="./stress_curves/freq_sweep/biaxial.out"
s="./stress_curves/freq_sweep/shear.out"
ss="./stress_curves/freq_sweep/sshear.out"

set origin 0.0,0.5
set xlabel "Frequency"
set ylabel "Storage Modulus"
set log y
set log x
plot m using 1:2 with lines title 'uniaxial', \
     p using 1:2 with lines title 'equibiaxial', \
     s using 1:2 with lines title 'shear', \
     ss using 1:2 with lines title 'simple shear' 
set origin 0.0,0.0
set xlabel "Frequency"
set ylabel "Loss Modulus"
set log y
set log x
plot m using 1:3 with lines title 'uniaxial', \
     p using 1:3 with lines title 'equibiaxial', \
     s using 1:3 with lines title 'shear', \
     ss using 1:3 with lines title 'simple shear' 
unset multiplot
#plot m using 1:3  with lines title 'passive',m using 1:4 with lines title 'active'
#plot    m using 2:3  with lines title 'Lissajous'
