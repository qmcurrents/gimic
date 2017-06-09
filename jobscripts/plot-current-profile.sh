#!/bin/bash

function createPlot() {

echo "Creating plot, please wait..."

gnuplot << EOF

# colour and style definitions:
#cyan: rgb 0,255,255 #00FFFF
#blue: rgb 30,70,255 #1E46FF
#green: rgb 0,178,0  #00B200

# diatropic (green)
set style line 1 lt 1 lw 5 lc rgb "#007F00" 
# paratropic (blue)
set style line 2 lt 3 lw 5 lc rgb "#1E46FF"
# vertical lines (cyan)
set style line 3 lt 1 lw 2 lc rgb "#00DCFF"
# vertical zero line
set style line 4 lt 1 lw 5 lc rgb "#000000" 

set format x "%5.2f"
set format y "%5.2f"
unset label
set xlabel "Distance [bohr]"
set ylabel "dJ/dx [nA/T / bohr]"

set key $plotkey 

# offset from the geometrical centre
offset=${offsets[0]}
#offset=-${offsets[0]}  # <-- for the weird bond 23.9 in BNBN
max_x=$max_x
set xrange [-offset:max_x-offset]

set terminal postscript eps enhanced color size 4.90, 2.70 'Helvetica' 22

# draw a vertical line at x=0
set arrow from 0, graph 0 to 0, graph 1 nohead ls 4

# Pseudo loop for the case of x=0 at the bond
# calculate how far the centre of the ring is from the "in" value and offset it by the default offset (of the bond)
if ( ${offsets[1]}!=0) set arrow from (${offsets[1]}), graph 0 to (${offsets[1]}), graph 1 nohead ls 4
if ( ${offsets[2]}!=0) set arrow from (${offsets[2]}), graph 0 to (${offsets[2]}), graph 1 nohead ls 4

# draw a vertical line at x=xin
#xin=offset;
#set arrow from xin,graph 0 to xin,graph 1 nohead ls 4
#set arrow from 3, graph 0 to 3, graph 1 nohead
#xin=1.5-offset
#set arrow from xin,graph(0,0) to xin,graph(1,1) nohead ls 3
#xin=2.54-offset
#set arrow from xin,graph(0,0) to xin,graph(1,1) nohead ls 3
#xin=4.22-offset
#set arrow from xin,graph(0,0) to xin,graph(1,1) nohead ls 3

#set output "current-profile.eps"
#plot "current_profile.dat" u ($1-offset):2 w l lw 5 notitle

set output "current_profile_$atom1.$atom2-dia-para.eps"
plot "current_profile.dat" u (\$1-${offsets[0]}):3 w l ls 1 title "Diatropic", "current_profile.dat" u (\$1-${offsets[0]}):4 w l ls 2 title "Paratropic"
# for the weird bond 23.9 in BNBN: 
#plot "current_profile.dat" u (\$1+${offsets[0]}):3 w l ls 1 title "Diatropic", "current_profile.dat" u (\$1+${offsets[0]}):4 w l ls 2 title "Paratropic"

EOF

display current_profile_$atom1.$atom2-dia-para.eps &
}

function centroid() { 
        atoms="$@"
        for i in $atoms 
        do
            awk -v i=$i '{  
                        if (NR == i+1) print $0 
                        }' coord 
        done | awk '{  
                     x+=$1; 
                     y+=$2; 
                     z+=$3; 
                     n++
                  } 
                  END { 
                     print x/n, y/n, z/n ; 
                 }';
}

function startCentroid() {
echo "DEFINE THE POINT x=0 ON THE PLOT"
echo "Enter the indices of the atoms according to the coord file and then press ENTER"
read atoms

echo "centroid:"
cent=$(centroid $atoms)
echo $cent

Cx=$( echo $cent | awk '{print $1} ' )
Cy=$( echo $cent | awk '{print $2} ' )
Cz=$( echo $cent | awk '{print $3} ' );

}

function bondCentroid() {
wrkdir="$1"
if [ -z "$wrkdir" ]
then
    wrkdir=$(pwd)
fi
# Find the name of the work directory
#dir=$(pwd) # take the whole path
dirname=${wrkdir##*/} # pick the name of the current directory only

echo "Dirname: $dirname"

atom1=$(echo $dirname | awk -F'[^0-9]*' '{print $2}') # pick the first number from the name of the directory 
atom2=$(echo $dirname | awk -F'[^0-9]*' '{print $3}') # pick the second number from the name of the directory

bond_cent=$(centroid $atom1 $atom2)

# Coordinates of the centre:
BOND=$(centroid $atom1 $atom2)
Bx=$( echo $BOND | awk '{print $1} ' )
By=$( echo $BOND | awk '{print $2} ' )
Bz=$( echo $BOND | awk '{print $3} ' )

echo "Bond: $atom1 -- $atom2"
echo "Bond centre: ($Bx, $By, $Bz)"

start=$( sed -n -e 's/^.*in=//p' calculation.dat | awk '{print $1}')
echo "in = $start"

}

function calculateOffset() {
# $distance is the distance between the beginning of the int plane and the bond - it should not be changed any further by the script
distance=$( awk -v Cx=$Cx -v Cy=$Cy -v Cz=$Cz -v Bx=$Bx -v By=$By -v Bz=$Bz BEGIN'{
            dist=sqrt( (Bx-Cx)*(Bx-Cx) + (By-Cy)*(By-Cy) + (Bz-Cz)*(Bz-Cz) ); print dist}' )
echo "Distance centre to bond: $distance"

offset=$( awk -v dist=$distance -v start=$start 'BEGIN{print start-dist}' )
}

function calculateDistanceAfter() {
distanceAfter=$( awk -v Cx=$Cx -v Cy=$Cy -v Cz=$Cz -v Bx=$Bx -v By=$By -v Bz=$Bz BEGIN'{
                  dist=sqrt( (Bx-Cx)*(Bx-Cx) + (By-Cy)*(By-Cy) + (Bz-Cz)*(Bz-Cz) ); print dist}' )
echo "DistanceA centre to bond: $distanceAfter"
offsetAfter=$( awk -v dist=$distance -v offset=$offset -v after=$distanceAfter 'BEGIN{print dist+after}' )
}


function calculateDistanceBefore() {
distanceBefore=$( awk -v Cx=$Cx -v Cy=$Cy -v Cz=$Cz -v Bx=$Bx -v By=$By -v Bz=$Bz BEGIN'{
                  dist=sqrt( (Bx-Cx)*(Bx-Cx) + (By-Cy)*(By-Cy) + (Bz-Cz)*(Bz-Cz) ); print dist}' )
echo "DistanceB centre to bond: $distanceBefore"
offsetBefore=$( awk -v dist=$distance -v start=$start 'BEGIN{print dist-start}' )
#offsetBefore=$( awk -v distB=$distanceBefore 'BEGIN{ print -distB}' )
}




offsets=( )
offsets[0]=0
offsets[1]=0
offsets[2]=0

echo
startCentroid;
echo $atoms
echo "Vertical axis atoms=($atoms)" >> calculation.dat

bondCentroid;
calculateOffset;
echo "Cx=$Cx Cy=$Cy Cz=$Cz" >> calculation.dat
echo "Bx=$Bx By=$By Bz=$Bz" >> calculation.dat

idx=0;
offsets[$idx]=$offset  # save the original value of the offset because it will be modified
echo "offsets[$idx] = ${offsets[$idx]}";
echo "offset=${offsets[$idx]}" >> calculation.dat

max_x=$(tail -n 1 current_profile.dat | awk 'END {print $1}')
printf  "\nWould you like to define other vertical lines?\nPress [S] if the line lies at x < 0 bohr.\nPress [L] if the line lies at x > 0 bohr.\n Or press [ENTER] to skip.\n"
read accept;

while [ ! -z $accept ]
do
    startCentroid;
    idx=$(($idx+1));
    if [ "$accept" == "L" ]
    then 
	# for x = 0 at a bond between two rings
    calculateDistanceAfter;
	offsets[$idx]=$offsetAfter
    elif [ "$accept" == "S" ]
    then
    calculateDistanceBefore;
	# for x = 0 at a bond between two rings
	offsets[$idx]=$offsetBefore
    fi
    echo "offsets[$idx] = ${offsets[$idx]}";
printf  "\nWould you like to define other vertical lines?\nPress [S] if the line lies at x < 0 bohr.\nPress [L] if the line lies at x > 0 bohr.\n Or press [ENTER] to skip.\n"
    read accept;
done
# NOTE: the first plane will be before the bond, and the second after the bond -> NEEDS TO BE FIXED!

plotkey=default
createPlot;

printf  "\nDo you accept the $plotkey position of the key in the plot?\nPress [A] to move it above the plot\nPress [B] to move it at the bottom right corner of the plot\nOr press [ENTER] to exit.\n"
read accept;
if [ ! -z $accept ] && [ $accept == "A" ]
then
    plotkey="above"
    createPlot;
fi
if [ ! -z $accept ] && [ $accept == "B" ]
then
    plotkey="right bottom"
    createPlot;
fi


