set xlabel " time"
set ylabel " gamma"
m="./simple_shear.out"
#set terminal qt
set terminal png
set output 'simple_shear_rheology.png'
#set terminal postscript eps enhanced color font 'Helvetica,12'
#set output 'passive_vs_active_shear.eps'
#set nokey
set key bottom right
set grid
#set log y
#set log x
#plot m using 1:2 wit lines title 'gamma'
#plot m using 1:2  with lines title 'gamma'
#plot    m using 1:3  with lines title 'Sxy'
plot m using 1:2 with lines title 'storage modulus', m using 1:3 with lines title 'loss modulus' 
#plot m using 1:3  with lines title 'passive',m using 1:4 with lines title 'active'
#plot    m using 2:3  with lines title 'Lissajous'
