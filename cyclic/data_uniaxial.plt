set xlabel " stretch"
set ylabel " stress"
m="./uniaxial.out"
set terminal postscript eps enhanced color font 'Helvetica,12'
set output 'uniaxial.eps'
#set nokey
set key bottom right
set grid
#set log y
#set log x
plot m using 2:3 wit lines title 'uniaxial'
#plot m using 1:3  with lines title 'passive',m using 1:4 with lines title 'active'
#plot    m using 2:3  with lines title 'Lissajous'
