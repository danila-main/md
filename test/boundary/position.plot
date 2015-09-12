set term X11
#set output "profit.png"

set multiplot layout 2, 1 title "Two noninteracting particles move along X: T=10, dt=1e-11"

set xlabel 'Step no' offset 0,0.9
set ylabel 'Position' offset 3.5,0
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



set xlabel 'Step no' offset 0,0.9
set ylabel 'Distance' offset 3.5,0
set yrange [0:0.55]
#set key right top
set nokey
plot 'part_pos' u 1:($4/0.00170999) with lines lt 1 lw 1.5 notitle,\
	  (0.5*0.00170999/0.00170999) with lines lt 1 lc rgb"black" lw 1 
	  #t '1/2 Cell edge'


pause -1
