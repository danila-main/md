#set term wxt
#set output "profit.png"
set terminal pngcairo size 1024,768 crop enhanced font 'Verdana,16'
set output "move_b_b_pos.png"

set xlabel 'x/L, 10^{2}' 
set ylabel 'y/L, 10^{2}'
set size square
set yrange [49.82:50.02]
set mxtics
set grid
set grid mxtics
set nokey
set title 'H = 5x10^{4} G, {/Symbol D}t = 10^{-11}'
plot 'part_pos_b' u ($2*100):($3*100) w p lt 7 ps 0.2 lc rgb"red"  

#pause -1
