#!/bin/bash


# Function definitions:
source  /appl/gimic-fork/gimic/jobscripts/functions-def

function critical_dia() { 
awk 'BEGIN { 
      prev=10; 
      prev2=0; 
  }
  { 
      diff1 = prev-prev2;
      distance2=$1
      if ( diff1 < 1e-6 ) { 
          diff2 = $3 - prev;
          if ( diff2 > 1e-6) {
              printf("%.2f\n",distance); 
          }
          if ((prev > 1e-6)&&($3 < 1e-6)) {
              printf("%.2f\n",distance2);
          }
      } 
      prev2=prev; 
      prev=$3; 
      distance=$1
    }' $wrkdir/current_profile.dat
};

function critical_para() {
awk 'BEGIN{ 
      prev=-10;  
      prev2=0; 
  }
  {
      diff1 = prev-prev2;
      distance2=$1
      if ( diff1 > 1e-6 ) { 
          diff2 = $4 - prev;
          if ( diff2 < 1e-6) {
              printf("%.2f\n",distance); 
          }
      }
      if ((prev > -1e-6)&&($4 < -1e-6)) {
              printf("%.2f\n",distance2); 
          }

      prev2=prev; 
      prev=$4; 
      distance=$1
    }' $wrkdir/current_profile.dat

}; 

function critical_net() {
awk 'BEGIN{ 
      prev=$2;  
  }
  {
      if ( (( prev < 0 )&&($2 > 0 )) || (( prev > 0 )&&($2 < 0 )) ) { 
              printf("%.2f\n",$1); 
          }
      prev=$2; 
    }' $wrkdir/current_profile.dat
}; 

function coords() {
# usage: 
# echo $x0 $x1 | coords $length 
# note: $x0 and $x1 contain all three coordinates of the points
#
#  We need to calculete the length of a small piece parallel to the (x,y,0) plane
#  It is calculated using the angle beta between this parallel line and the whole line, which is tilted in the general case
#	beta = arcsin(dz / totLength); 
#  Then the length of the small piece is the cosine of this angle beta times the length of the piece of the whole line:
#	dist = cos(beta)*lngt;
#  Using cos(arcsin(x)) = sqrt(1 - x^2):
#  dist = sqrt(1 - dz*dz/totLength/totLength)*lngt

awk -v lngt=$1 -v totLength=$length '{
        # function arcsin(x) { return atan2(x,sqrt(1-x*x)) } 
	split($0,crd); 
	dx = (crd[1]-crd[4]); 
	dy = (crd[2]-crd[5]); 
	dz = (crd[3]-crd[6]); 
	alpha = atan2(dy, dx);
#print ("Total length = ", totLength);
#print ("Local length = ", lngt);
#print ("dz =", dz);
#print ("1 - dz*dz/totLength/totLength = ", 1 - dz*dz/totLength/totLength);

	dist = sqrt(1 - dz*dz/totLength/totLength)*lngt;

#print "\n dist = ", dist;
#print "cos(alpha) = ", cos(alpha);
#print "sin(alpha) = ", sin(alpha);
	xcoord = -cos(alpha)*dist + crd[1];
	ycoord = -sin(alpha)*dist + crd[2];
	print "\tCoords: ( " xcoord, ";", ycoord " )"; 
    }'
}

# When writing to the output file current_profile_#.#.txt this is in use: 
function coordsOut() 
{
#    awk -v alpha=$alpha -v x0x=$x0x -v x0y=$x0y -v lngt=$1 -v dz=$dz'{
    awk -v alpha=$alpha -v x0x=$x0x -v x0y=$x0y -v dz=$dz -v totLength=$length '{
	lngt=$1;
#	dist = sqrt(1 - dz*dz/totLength/totLength)*lngt;
        dist=lngt
	xcoord = -cos(alpha)*dist + x0x;
	ycoord = -sin(alpha)*dist + x0y;

#	printf("% 9.6f  % 9.6f  % 9.6f  % 9.6f  % 9.6f\n", xcoord, ycoord, $2+dz, $3+dz, $4+dz); 
	printf("% 9.6f  % 9.6f  % 9.6f  % 9.6f  % 9.6f\n", xcoord, ycoord, $2, $3, $4); 
    }' $wrkdir/current_profile.dat
};


###################################################################################################################

wrkdir="$1"
if [ -z "$wrkdir" ]
then
    wrkdir=$(pwd)
fi
checkIfEmpty wrkdir $wrkdir

# Find the name of the work directory
#dir=$(pwd) # take the whole path
dirname=${wrkdir##*/} # pick the name of the current directory only
checkIfEmpty dirname $dirname
#echo Working directory: $wrkdir
#echo Dirname: $dirname

# create empty files for the output
cat /dev/null > $wrkdir/profile-points.out
cat /dev/null > $wrkdir/profile-points-paraview.dat
#cat /dev/null > $wrkdir/profile-points-paraview.py

echo "Calculating the coordinates of the starting point of integration"
echo

(cd $wrkdir && gimic --dryrun gimic.0.inp > /dev/null ) 

x0A=$(awk '/X/{i++}i==1{print $2, $3, $4; exit}' "$wrkdir"/grid.xyz ) 

echo "Calculating the coordinates of the end point of integration"
last=$(( $(ls "$wrkdir"/gimic*inp | wc -l) -1 )) &&

( cd $wrkdir && gimic --dryrun gimic.$last.inp > /dev/null )

x1A=$(awk '/X/{i++}i==1{print $2, $3, $4; exit}' < "$wrkdir"/grid.xyz) 

# check which point lies higher to determine the start and end
crd0=(${x0A// / })
crd1=(${x1A// / })
x0x=$(A2bohr ${crd0[0]})
x0y=$(A2bohr ${crd0[1]})
x0z=$(A2bohr ${crd0[2]})
x1z=$(A2bohr ${crd1[2]})
dz=$(awk -v z0=$x0z -v z1=$x1z 'BEGIN{dz = z0-z1; print dz;}')

echo dz = $dz

length=$( echo $x0A $x1A | point_dist )
length=$(A2bohr $length)
echo "Length of the integration plane: "  $length, "bohr" >> $wrkdir/profile-points.out

#echo
#echo x0 [Å] = $x0A
#echo x1 [Å] = $x1A

x0=$(A2bohr $x0A)
x1=$(A2bohr $x1A)
echo
#echo x0 [bohr] = $x0
#echo x1 [bohr] = $x1

printf "\nin"  >> $wrkdir/profile-points.out
echo $x0 $x1 | coords 0  >> $wrkdir/profile-points.out
printf "out"  >> $wrkdir/profile-points.out
echo $x0 $x1 | coords $length  >> $wrkdir/profile-points.out

printf "\nDiatropic critical points:\n"  >> $wrkdir/profile-points.out
dia=$(critical_dia)
for value in $dia
do
    printf $value >> $wrkdir/profile-points.out
    echo $x0 $x1 | coords $value >> $wrkdir/profile-points.out
#    awk -v value=$value '{if ($1 == value) {print $3}}'  $wrkdir/current_profile.dat
done

#echo "DIA DONE"

printf "\nParatropic critical points\n"  >> $wrkdir/profile-points.out
para=$(critical_para)
for value in $para
do
    printf $value >> $wrkdir/profile-points.out
    echo $x0 $x1 | coords $value >> $wrkdir/profile-points.out
#    awk -v value=$value '{if ($1 == value) {print $4}}'  $wrkdir/current_profile.dat
done

#echo "PARA DONE"

printf "\nThe net current changes sign at the points:\n" >> $wrkdir/profile-points.out
net=$(critical_net)
for value in $net 
do
    printf $value >> $wrkdir/profile-points.out
    echo $x0 $x1 | coords $value >> $wrkdir/profile-points.out
 #   awk -v value=$value '{if ($1 == value) {print $2}}'  $wrkdir/current_profile.dat
done

cat $wrkdir/profile-points.out
echo


printf "\nOUTPUT RELATED TO PARAVIEW\n\n"
echo "Origin of the clipping plane:"

# Calculate the centre of the bond
# atom1=$( echo $dirname | awk -F'[^0-9]*' '$0=$2' ) # Same as the method below
# atom2=$( echo $dirname | awk -F'[^0-9]*' '$0=$3' )
atom1=$(echo $dirname | awk -F'[^0-9]*' '{print $2}') # pick the first number from the name of the directory 
atom2=$(echo $dirname | awk -F'[^0-9]*' '{print $3}') # pick the second number from the name of the directory

centroid $atom1 $atom2

printf "\n# Clipping plane origin:\n# " >> $wrkdir/profile-points-paraview.dat
echo "(" $( centroid $atom1 $atom2 ) ")" >> $wrkdir/profile-points-paraview.dat

# find the normal to the clipping plane in Paraview:
echo; echo "Normal to the clipping plane:"
echo $x0 $x1 | normalPlane 1
echo

printf "\n# Normal to the clipping plane:\n# " >> $wrkdir/profile-points-paraview.dat
echo $x0 $x1 | normalPlane 1 >> $wrkdir/profile-points-paraview.dat
echo >> $wrkdir/profile-points-paraview.dat

# Calculate the zero points in the Gnuplot profile

alpha=$( (echo  $x0 $x1 | angle) ) # calls the angle() function to find the angle of the slope of the line
#echo "alpha = " $alpha

delta=$( sed -n -e 's/^.*delta=//p' $wrkdir/calculation.dat | awk '{print $1}')
checkIfEmpty delta $delta
#echo "Slice thickness: " $delta
#echo

#echo x0 x,y: $x0x $x0y

# Write the coordinates and currents data
echo " X  Y  T  D  P" > $wrkdir/$dirname.txt 
echo "written to file  $wrkdir/$dirname.txt"
#echo $x0 $x1 | coordsOut $length  >> $wrkdir/$dirname.txt
coordsOut >> $wrkdir/$dirname.txt

## Write the Python functions for the Paraview visualizations
## -> the line along which the plane cuts the molecule and the zero points
#awk '{  
#        if ($1 == "in" ) {
#            pt="point_" $1; 
#            print pt " = Sphere()"; 
#            print pt ".Center = [" $4 ", " $6 ", 0.0]"; 
#            print pt ".Radius = 0.08"; 
#            print pt ".ThetaResolution = 20"; 
#            print pt ".PhiResolution = 20"; 
#            print "Show(" pt ")\n"
#
#            print "plane = Line()"; 
#            print "plane.Point1 = [" $4 ", " $6 ", 0.0]"; 
#        } 
#    }' profile-points.out >> $wrkdir/profile-points-paraview.py
#
#awk '{  
#        if ($1 == "out" ) { 
#            print "plane.Point2 = [" $4 ", " $6 ", 0.0]";  
#            print "Show(plane)\n"
#            
#
#            pt="point_" $1; 
#            print pt " = Sphere()"; 
#            print pt ".Center = [" $4 ", " $6 ", 0.0]"; 
#            print pt ".Radius = 0.08"; 
#            print pt ".ThetaResolution = 20"; 
#            print pt ".PhiResolution = 20"; 
#            print "Show(" pt ")\n"
#        } 
#    }' profile-points.out  >> $wrkdir/profile-points-paraview.py
#
#awk '{  
#        if ($1 ~ /^[0-9]/ ) {
#            idx++;
#            pt="point_" idx; 
#            print pt " = Sphere()"; 
#            print pt ".Center = [" $4 ", " $6 ", 0.0]"; 
#            print pt ".Radius = 0.08"; 
#            print pt ".ThetaResolution = 20"; 
#            print pt ".PhiResolution = 20"; 
#            print "Show(" pt ")\n"
#        } 
#    }' profile-points.out >> $wrkdir/profile-points-paraview.dat
#
