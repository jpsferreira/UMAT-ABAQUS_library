set xlabel "Generation"
set ylabel "Fitness"
m="./plot_ga.out"
#set terminal x11 0
set terminal postscript eps enhanced color
set output 'ga_soft_tissue.eps'
#set nokey
set key top left
set grid
set title 'Soft tissue: Genetic algorithm fitness'
plot m using 1:2  with lines title 'Avg', m using 1:3  with lines title 'Best'