set format x "%4.1f"
set format y "%4.1f"
set format z "%4.1f"
unset label
set xlabel "Distance [bohr]"
set ylabel "dJ/dx [nA/T / bohr]"

#set key above

set terminal postscript eps enhanced color 'Helvetica' 22

set xrange [-0.5:9]

# draw a vertical line at x=in
#set arrow from 2.4,graph(0,0) to 2.4,graph(1,1) nohead lw 3
set arrow from 0,graph(0,0) to 0,graph(1,1) nohead lw 3

offset=0.6

set output "current-profile.eps"
plot "current_profile.dat" u ($1-offset):2 w l lw 5 notitle
# plot "current_profile.dat" u 1:2 w l lw 5 notitle

set output "current-dia-para.eps"
plot "current_profile.dat" u ($1-offset):3 w l lc 3 lw 5 title "Diatropic", "current_profile.dat" u ($1-offset):4 w l lc 1 lw 5 title "Paratropic"
# plot "current_profile.dat" u 1:3 w l lc 3 lw 5 title "Diatropic", "current_profile.dat" u 1:4 w l lc 1 lw 5 title "Paratropic

