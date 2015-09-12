#set term wxt
#set output "profit.png"
set terminal pngcairo size 1024,768 crop enhanced font 'Verdana,16'
set output "move_b_en.png"

set xlabel 'Step no' 
set ylabel '1-E/E_{0}'
set mxtics
set grid
set grid mxtics
#set yrange [0:0.55]
#unset yrange
#unset size
set logscale y
set format y "10^{%L}"
set key right top
#set nokey
#plot 'part_pos' u 1:($4/0.00170999) with lines lt 1 lw 1.5 notitle,\
#	  (0.5*0.00170999/0.00170999) with lines lt 1 lc rgb"black" lw 1 
	  #t '1/2 Cell edge'
plot 'energy_vv' u 1:(abs($4)) with lines lt 1 lw 1.5 lc rgb"blue" t 'Verlet', 'energy_b' u 1:(abs($4)) with lines lt 1 lw 1.5 lc rgb"red" t 'Current'

#pause -1
