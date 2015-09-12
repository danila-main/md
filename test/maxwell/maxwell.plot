#
set term wxt
set output "profit.png"
set multiplot layout 2, 1 title "Number of particles: 1000, T=10"

bin_width = 0.1e6; # edit this 
bin_number(x) = floor(x/bin_width)
rounded(x) = bin_width * ( bin_number(x) + 0.5 )
UNITY = 1
# column number of data to be histogrammed is here assumed to be 1
# - change $1 to another column if desired
#plot 'el_V_dist' u (rounded($1)):(UNITY) t 'data' smooth frequency w histeps
set xlabel 'Velocity' offset 0,0.9
set ylabel 'Number of particles' offset 3.5,0
set boxwidth 0.7 relative
set style fill solid 0.5 noborder
plot 'exp_dist' u (rounded($1)):(UNITY) smooth frequency w boxes t 'Generated distribution Vx', \
	 'th_dist' u 1:2 with points ps 0.7 pt 7 t 'Maxwellian distribution Vx'
#	 'exp_dist_adj' u (rounded($1)):(UNITY) smooth frequency w boxes fs fill empty border lt 1 lc rgb"blue" \
#					  t 'Distribuiton of adjusted Vx'



#set title "Vector distribution"
plot 'exp_dist' u (rounded($4)):(UNITY) smooth frequency w boxes t 'Generated distribution V', \
	 'th_dist' u 3:4 with points ps 0.7 pt 7 t 'Maxwellian distribution V'
#	 'exp_dist_adj' u (rounded($4)):(UNITY) smooth frequency w boxes fs fill empty border lt 1 lc rgb"blue" \
#					  t 'Distribuiton of adjusted V'


pause -1
