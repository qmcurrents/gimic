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

function startCentroid() {
Cx=$( sed -n -e 's/^.*Cx=*//p' calculation.dat | awk '{print $1}')
Cy=$( sed -n -e 's/^.*Cy=*//p' calculation.dat | awk '{print $1}')
Cz=$( sed -n -e 's/^.*Cz=*//p' calculation.dat | awk '{print $1}')

if [ ! -z "$Cx" ] && [ ! -z "$Cy" ] && [ ! -z "$Cz" ]
then
    echo "Start centroid: ($Cx; $Cy; $Cz)"
else

    echo "Enter the indices of the atoms according to the coord file"
    read atoms
    echo $atoms

    echo "Centroid:"
    cent=$(centroid $atoms)
    echo $cent

    Cx=$( echo $cent | awk '{print $1} ' )
    Cy=$( echo $cent | awk '{print $2} ' )
    Cz=$( echo $cent | awk '{print $3} ' );

    # echo "Cx=$Cx Cy=$Cy Cz=$Cz" >> calculation.dat
fi

}

function bondCentroid() {
Bx=$( sed -n -e 's/^.*Bx=*//p' calculation.dat | awk '{print $1}')
By=$( sed -n -e 's/^.*By=*//p' calculation.dat | awk '{print $1}')
Bz=$( sed -n -e 's/^.*Bz=*//p' calculation.dat | awk '{print $1}')

if [ ! -z "$Bx" ] && [ ! -z "$By" ] && [ ! -z "$Bz" ]
then
    echo "Bond centroid: ($Bx; $By; $Bz)"
else
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

    echo "Bond: $atom1 $atom2"
    echo "Bond centre: ($Bx, $By, $Bz)"
    #echo "Bx=$Bx By=$By Bz=$Bz" >> calculation.dat
fi
}

function calculateDistOff() {
offset=$( sed -n -e 's/^.*offset=*//p' calculation.dat | awk '{print $1}')

if [ ! -z "$offset" ]
then
    echo "Offset = $offset"
else
    start=$( sed -n -e 's/^.*in=//p' gimic.0.inp | awk '{print $1}')
    echo "in = $start"

    distance=$( awk -v Cx=$Cx -v Cy=$Cy -v Cz=$Cz -v Bx=$Bx -v By=$By -v Bz=$Bz BEGIN'{dist=sqrt( (Bx-Cx)*(Bx-Cx) + (By-Cy)*(By-Cy) + (Bz-Cz)*(Bz-Cz) ); print dist}' )
    echo "Distance centre to bond: $distance"

    offset=$( awk -v dist=$distance -v start=$start 'BEGIN{print start-dist}')
    echo "Offset: $offset";
fi
}

startCentroid
bondCentroid
calculateDistOff

intervals=( );
idx=0;
interval="interval"
lower=0;
upper=0;

while [ ! -z "$interval" ] 
do
    echo "Enter interval boundaries and then press Enter:"
    read interval
    intervals[$idx]=$( echo $interval )
    idx=$(($idx+1));
done

idx=$(($idx-1));

# Find the starting value and set it as an offset:


for (( i=0; i<$idx; i++ ))
do
    lower=$(echo ${intervals[$i]} | awk '{print $1}');
    upper=$(echo ${intervals[$i]} | awk '{print $2}');
    awk -v lower=$lower -v upper=$upper -v offset=$offset '{  
       if (($1 >= lower) && ($1 <= upper)) { 
           total+=$2; 
           dia+=$3; 
           para+=$4; } 
       } 
       END { 
           printf("%.2f & %.2f & %.2f & %.2f & %.2f \\\\ \n", lower-offset, upper-offset, total, dia, para); 
       } ' current_profile.dat ; 
done

