set xlabel " Frequency"
set ylabel " G2"
m="./data.out"
n="./exp_data.out"
set terminal x11 0
#set terminal png
#set output 'passive_vs_active_shear.png'
#set terminal postscript eps enhanced color font 'Helvetica,12'
#set output 'passive_vs_active_shear.eps'
#set nokey
set key bottom right
set grid
set log y
set log x
plot m using 1:2  with lines title 'G', m using 1:3  with lines title 'tan', \
     n using 1:2  title 'Gexp', n using 3:4 title 'tanexp'
