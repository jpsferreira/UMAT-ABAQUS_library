set xlabel "Amount of Stretch"
set ylabel "PK1 Stress (Nominal Stress)"
m="./plot_.out"
r="./plot_r2.out"
rr=system("head -n1 plot_r2.out | awk '{print $1}'")
p1=system("head -n2 plot_r2.out | tail -n1 | awk '{print $1}'")
p2=system("head -n2 plot_r2.out | tail -n1 | awk '{print $2}'")
p3=system("head -n2 plot_r2.out | tail -n1 | awk '{print $3}'")
p4=system("head -n2 plot_r2.out | tail -n1 | awk '{print $4}'")
#set terminal x11 0
set terminal postscript eps enhanced color
set output 'soft_tissue.eps'
#set nokey
set key top left
set grid
set label 1 gprintf("R^2 = %3.4f",rr) right at graph 0.2, graph 0.8
set label 2 gprintf("C_{10} = %3.4E",p1) right at graph 0.2, graph 0.7
set label 3 gprintf("k_{1} = %3.4E",p2) right at graph 0.2, graph 0.65
set label 4 gprintf("k_{2} = %3.4E",p3) right at graph 0.2, graph 0.6
set label 5 gprintf("k_{disp} = %3.4E",p4) right at graph 0.2, graph 0.55
set title 'Soft tissue: Uniaxial tension fitting'
plot m using 1:2 title 'Exp', m using 1:3  with lines title 'Fit sxx'
