#set term wxt
#set output "profit.png"
set terminal pngcairo size 1024,768 enhanced font 'Verdana,16'
set output "move_e.png"
set multiplot layout 2, 1 title "Two interacting particles move along X: T=0, dt=1e-13"

set xlabel 'Step no' 
set ylabel 'Position'
set yrange [0:1.1]
set mxtics
set grid
set grid mxtics
set key right bottom
plot 'part_pos' u 1:($2/0.00170999) with lines lt 1 lc rgb"red" lw 1.5 t 'Electron', \
	 'part_pos' u 1:($3/0.00170999) with lines lt 1 lc rgb"blue" lw 1.5 t 'Proton', \
      (1.0*0.00170999/0.00170999) with lines lt 1 lc rgb"black" lw 1 notitle , \
	  (0.5*0.00170999/0.00170999) with lines lt 1 lc rgb"black" lw 1 notitle 
	  #0.00170999 with lines lt 1 lc rgb"black" lw 1 t 'Cell edge', \
	  #3.0/4.0 * 0.00170999 with lines lt 1 lc rgb"black" lw 1 notitle



set xlabel 'Step no' 
set ylabel '1-E/E_{0}'
#set yrange [0:0.55]
unset yrange
set logscale y
set format y "10^{%L}"
#set key right top
set nokey
#plot 'part_pos' u 1:($4/0.00170999) with lines lt 1 lw 1.5 notitle,\
#	  (0.5*0.00170999/0.00170999) with lines lt 1 lc rgb"black" lw 1 
	  #t '1/2 Cell edge'
plot 'energy' u 1:4 with lines lt 1 lw 1.5 

#pause -1
