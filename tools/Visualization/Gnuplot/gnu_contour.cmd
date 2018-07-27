# usage: gnuplot gnu_contour.cmd
# purpose: get a quick visualization of jmod.txt after running cdens

set format x "%4.1f"
set format y "%4.1f"
set format z "%4.1f"
set contour base
#                    startvalue, increment, end
set cntrparam levels incremental 0.00, 0.001, 5.50
set style data lines 

unset surface
set view 0,0
unset key
set size square
set zrange[*:*]
set xrange[*:*]
set yrange[*:*]

# splot 'jmod.txt' using ($1):($3):($4)  
# splot 'jmod.txt' using ($2):($3):($4)  
  splot 'jmod.txt' using ($1):($2):($4)  

pause -1

# save output as eps file 
set terminal X11
# set terminal postscript eps enhanced monochrome 'Helvetica' 22
set terminal postscript eps enhanced color 'Helvetica' 22
set output 'contour.eps'

replot
pause -1

