#!/bin/bash

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

function A2bohr() { for val in $@; do awk -v a=$val 'BEGIN{ bohr=1.88971616463207; print a*bohr; }'; done; };
function bohr2A() { for val in $@; do awk -v a=$val 'BEGIN{ bohr=0.5291772490; print a*bohr; }'; done; };

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

function coords() {
awk -v dist=$1 '{
	split($0,crd); 
	dx=(crd[1]-crd[4]); 
	dy=(crd[2]-crd[5]); 
	dz=(crd[3]-crd[6]); 
	alpha=atan2(dy, dx);
	xcoord = -cos(alpha)*dist + $1;
	ycoord = -sin(alpha)*dist + $2;
	print "\tCoords: ( " xcoord, ";", ycoord " )"; 
    }'
}

# When writing to the output file current_profile_#.#.txt this is in use: 
function coordsOut() 
{
    awk -v alpha=$alpha -v x0x=$x0x -v x0y=$x0y '{
        dist=$1
	xcoord = -cos(alpha)*dist + x0x;
	ycoord = -sin(alpha)*dist + x0y;
	printf("% 9.6f  % 9.6f  % 9.6f  % 9.6f  % 9.6f\n", xcoord, ycoord, $2, $3, $4); 
    }' $wrkdir/current_profile.dat
};

function normalPlane() 
{ 
# usage: 
# echo $x0 $x1 1 | normalPlane
# where the coordinates are in bohr
    # take the coordinates of x0 and x1
    # take a third point C that lies at (x0x, x0y, 1)
    # define vector a: C -> x0; vector b: C -> x1;
    # calculate their cross product to find the normal vector to the plane at point C

    awk '{
          split($0,crd);

	  cx = crd[1]; 
	  cy = crd[2]; 
	  cz = crd[7]; 

	  ax = cx - crd[1]; 
	  ay = cy - crd[2];
	  az = cz - crd[3];

	  bx = cx - crd[4]; 
	  by = cy - crd[5];
	  bz = cz - crd[6];

	  nx = -ay*bz + az*by ; 
	  ny = -az*bx + ax*bz ;
	  nz = -ax*by + ay*by ;

	  nmod=sqrt(nx*nx + ny*ny + nz*nz); 
	  print nx/nmod, ny/nmod, nz/nmod; 
    }'
};


function point_dist() 
{ 
awk '{
     split($0,crd);
     dx=crd[1]-crd[4]; 
     dy=crd[2]-crd[5]; 
     dz=crd[3]-crd[6]; 
     dist=sqrt(dx*dx + dy*dy + dz*dz);
     print dist; }'
};

function Delta() {
        awk '{ if (NR == 2) { delta=$1; print delta } }' current_profile.dat
    };  

function angle() 
{
    awk '{
         split($0,crd); 
	 dx=crd[1]-crd[4]; 
	 dy=crd[2]-crd[5]; 
	 dz=crd[3]-crd[6]; 

	 alpha=atan2(dy, dx); print alpha;

     }'
 };

wrkdir="$1"
if [ -z "$wrkdir" ]
then
    wrkdir=$(pwd)
fi
# Find the name of the work directory
#dir=$(pwd) # take the whole path
dirname=${wrkdir##*/} # pick the name of the current directory only


echo Working directory: $wrkdir
echo Dirname: $dirname
echo
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
# x_comp=$(awk -v x0x=${crd0[0]} -v x1x=${crd1[0]}  'BEGIN{ if (x0x<x1x) {print "1"; }  else {print "0"; } }')
# y_comp=$(awk -v x0y=${crd0[1]} -v x1y=${crd1[1]}  'BEGIN{ if (x0y<x1y) {print "1"; }  else {print "0"; } }')

length=$( echo $x0A $x1A | point_dist )
length=$(A2bohr $length)
echo; echo "Length of the integration plane: "  $length "bohr"

echo
echo x0 [Å] = $x0A
echo x1 [Å] = $x1A

x0=$(A2bohr $x0A)
x1=$(A2bohr $x1A)
echo
echo x0 [bohr] = $x0
echo x1 [bohr] = $x1

cat /dev/null > $wrkdir/profile-points.out
cat /dev/null > $wrkdir/profile-points-paraview.dat
cat /dev/null > $wrkdir/profile-points-paraview.py

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
done

printf "\nParatropic critical points\n"  >> $wrkdir/profile-points.out
para=$(critical_para)
for value in $para
do
    printf $value >> $wrkdir/profile-points.out
    echo $x0 $x1 | coords $value >> $wrkdir/profile-points.out
done

cat profile-points.out
echo

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
echo "Normal to the clipping plane:"
echo $x0 $x1 | normalPlane 1
echo

printf "\n# Normal to the clipping plane:\n# " >> $wrkdir/profile-points-paraview.dat
echo $x0 $x1 | normalPlane 1 >> $wrkdir/profile-points-paraview.dat
echo >> $wrkdir/profile-points-paraview.dat

# Calculate the zero points in the Gnuplot profile

alpha=$( (echo  $x0 $x1 | angle) ) # calls the angle() function to find the angle of the slope of the line
#echo "alpha = " $alpha

delta=$(Delta)
#echo "Slice thickness: " $delta
#echo

#echo x0 x,y: $x0x $x0y

# Write the coordinates and currents data
echo " X  Y  T  D  P" > $wrkdir/$dirname.txt 
echo "written to file  $wrkdir/$dirname.txt"
coordsOut >> $wrkdir/$dirname.txt

# Write the Python functions for the Paraview visualizations
# -> the line along which the plane cuts the molecule and the zero points
awk '{  
        if ($1 == "in" ) {
            pt="point_" $1; 
            print pt " = Sphere()"; 
            print pt ".Center = [" $4 ", " $6 ", 0.0]"; 
            print pt ".Radius = 0.08"; 
            print pt ".ThetaResolution = 20"; 
            print pt ".PhiResolution = 20"; 
            print "Show(" pt ")\n"

            print "plane = Line()"; 
            print "plane.Point1 = [" $4 ", " $6 ", 0.0]"; 
        } 
    }' profile-points.out >> $wrkdir/profile-points-paraview.py

awk '{  
        if ($1 == "out" ) { 
            print "plane.Point2 = [" $4 ", " $6 ", 0.0]";  
            print "Show(plane)\n"
            

            pt="point_" $1; 
            print pt " = Sphere()"; 
            print pt ".Center = [" $4 ", " $6 ", 0.0]"; 
            print pt ".Radius = 0.08"; 
            print pt ".ThetaResolution = 20"; 
            print pt ".PhiResolution = 20"; 
            print "Show(" pt ")\n"
        } 
    }' profile-points.out  >> $wrkdir/profile-points-paraview.py

awk '{  
        if ($1 ~ /^[0-9]/ ) {
            idx++;
            pt="point_" idx; 
            print pt " = Sphere()"; 
            print pt ".Center = [" $4 ", " $6 ", 0.0]"; 
            print pt ".Radius = 0.08"; 
            print pt ".ThetaResolution = 20"; 
            print pt ".PhiResolution = 20"; 
            print "Show(" pt ")\n"
        } 
    }' profile-points.out >> $wrkdir/profile-points-paraview.dat



