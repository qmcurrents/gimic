set nokey
set size square
set xtics 2.0
set ylabel "Current strength profile (in nA/T/bohr)" 
set xlabel "Extent in bohr "
#

plot "c_profile.out" u ($1*0.1):($2/0.1) w l lt -1 lw 3 t "dJ(x)/dx"

pause -1

