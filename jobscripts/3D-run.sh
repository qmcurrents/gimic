#!/bin/bash

function calcVertices() {

echo "calculating vertices" 


 cat coord | sed -n '2,/^\\$/p' |  awk 'BEGIN{
  			 infinity=10000;
                         minX=infinity; 
                         minY=infinity; 
                         minZ=infinity;  
                         maxX=-infinity; 
                         maxY=-infinity; 
                         maxZ=-infinity
		         offset=8} 
                         /\$/{next} {
                                 if ($1 < minX) {minX = $1;}; 
                                 if ($2 < minY) {minY = $2;};
                                 if ($3 < minZ) {minZ = $3;}; 
                                 if ($1 > maxX) {maxX = $1;}; 
                                 if ($2 > maxY) {maxY = $2;};
                                 if ($3 > maxZ) {maxZ = $3;}; 
                                 }  
                         END{print minX-offset, minY-offset, minZ-0.5*offset, maxX-minX+2*offset, maxY-minY+2*offset, maxZ-minZ+offset}' 
		     }
                         #END{print minX-8, minY-8, minZ-4, maxX, maxY, maxZ}' coord



################################################################################################################################
# DEBUG the min and max:

cat coord | sed -n "2,/^\\$/p" |  awk 'BEGIN{
  			 infinity=10000;
                          minX=infinity; 
                          minY=infinity; 
                          minZ=infinity;  
                          maxX=-infinity; 
                          maxY=-infinity; 
                          maxZ=-infinity} 
                          /\$/{next} {
                                  if ($1 < minX) {minX = $1;}; 
                                  if ($2 < minY) {minY = $2;};
                                  if ($3 < minZ) {minZ = $3;}; 
                                  if ($1 > maxX) {maxX = $1;}; 
                                  if ($2 > maxY) {maxY = $2;};
                                  if ($3 > maxZ) {maxZ = $3;}; 
                                  }  
                          END{print "Min x, y, z: ", minX, minY, minZ; print "Max x, y, z: ", maxX, maxY, maxZ}' 
 


if [ -d 3D ] # Does it exists
then
        printf "\n\n*** Directory 3D already exists.\n"
        echo "Enter [y] to overwrite or any key to exit."; read accept
        if [ -n "$accept" ] && [ "$accept" == "y" ]  # if the variable is not empty or equal to "y", remove the dir
        then
            rm -rf 3D/* # delete it if the user does not want to keep the existing dir
        else
            exit 
        fi
else
    mkdir 3D
fi
#printInputFile

#cp /hdd/gimic/jobscripts/gimic.cdens 3D/gimic.inp

vertices=$( calcVertices )
#echo Vertices: $vertices

vert=(${vertices// / })

#echo MINX = ${vert[0]}
#echo MINY = ${vert[1]}
#echo MINZ = ${vert[2]}


#echo ${vert[3]}
#echo ${vert[4]}
#echo ${vert[5]}

echo "Magnetic field direction"
# default along the Z axis
MFx=0.0
MFy=0.0
MFz=-1.0

echo "Do you accept the default MF orientation along the Z axis?"
echo "Press [n] to calculate the direction automatically or [e] to enter manually"

read accept;
if [ -n "$accept" ] 
then
    if [ "$accept" == "n" ]
    then
	#Bcoords=$(maximise_projection grid.xyz | sed -e 's#{#_#g; s#}#_#g; s#,#_#g' | awk -F [_] '{ {print $2, $3, $4} }')
	echo debug:
	maximise_projection coord.xyz 
	maximise_projection coord.xyz > 3D/field.dat
	MFx=$( cat 3D/field.dat | sed -e 's#{#_#g; s#}#_#g; s#,#_#g' | awk -F [_] '{ {print -$2} }')
	MFy=$( cat 3D/field.dat | sed -e 's#{#_#g; s#}#_#g; s#,#_#g' | awk -F [_] '{ {print -$3} }')
	MFz=$( cat 3D/field.dat | sed -e 's#{#_#g; s#}#_#g; s#,#_#g' | awk -F [_] '{ {print -$4} }')
    elif [ "$accept" == "e" ]
    then
	print "MFx = "; read MFx
	print "MFy = "; read MFy
	print "MFz = "; read MFz
    fi
fi


echo "Magnetic field vector coordinates: ($MFx; $MFy; $MFz)"

sedstring="s/@MINX@/${vert[0]}/; s/@MINY@/${vert[1]}/; s/@MINZ@/${vert[2]}/; s/@LX@/${vert[3]}/; s/@LY@/${vert[4]}/; s/@LZ@/${vert[5]}/; s/@MFx@/$MFx/; s/@MFy@/$MFy/; s/@MFz@/$MFz/;"

sed "$sedstring" /hdd/gimic/jobscripts/gimic.cdens > 3D/gimic.inp

#x0A=$(awk '/X/{i++}i==1{print $2, $3, $4; exit}' "$wrkdir"/grid.xyz ) 
#crd0=(${x0A// / })
#x0x=$(A2bohr ${crd0[0]})

