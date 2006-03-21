set contour base
#set hidden3d
set cntrparam bspline
set cntrparam levels incremental -2,0.1,2
#set cntrparam levels 15
#set noclabel
#show contour

#set term postscript enhanced color
#set out 'foo.ps'

set nokey 
set zrange [0:5]
splot 'JMOD.1.txt' not w l
pause -1

#set nosurface
#set key 3.7, 4
set label 'c' at 0,0
set label 'o' at -2.16,0
set label 'o' at 2.16,0
#set view 0,0,1.7
#set size square 1
#splot 'DIVJPLT.1.txt' using 3:1:4 not w l
#pause -1

#set size square 1
#set xrange [-4:4]
#set yrange [-4:4]
plot 'JVEC.1.txt' using 3:1:6:4 not w vector
pause -1
#plot 'NJVEC'  not w vector
#pause -1

#set size square 1
#set term table
#set output 'table'
#splot 'ampere' not w l
#set term postscript enhanced
#set output 'contour.ps'
#plot 'table' index 0:100 not w l
#pause -1
