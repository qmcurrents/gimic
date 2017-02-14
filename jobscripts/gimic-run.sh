#!/bin/bash

function centroid() { atoms="$@"; for i in $atoms; do awk -v i=$i '{ if (NR == i+1) print $0 }' coord; done | awk '{ x+=$1; y+=$2; z+=$3; n++} END{ print x/n, y/n, z/n; }'; }
function dista() { awk -v a1="$1" -v a2="$2" '{ if (NR == a1+1) { x+=$1;y+=$2;z+=$3;}; if (NR == a2+1) {x-=$1;y-=$2;z-=$3} } END{ dist=sqrt(x*x+y*y+z*z); print dist*0.5; }' coord ; }

echo
echo "Bond between the atoms:"
read atom1; read atom2
bond=$atom1","$atom2

echo "Enter the suffix for the directory name"
read suffix
suffix=$(printf _$suffix)
dirname=bond_int_$atom1.$atom2$suffix
echo $dirname

# If we are not in the current_profile folder, check if it exists. If it did, just go there. If it did not, create it and copy the necessary files
dir=$(pwd | grep "bond_int")
#echo $dir
if [ -z $dir ] # if the path does not contain "current_profile", $dir is empty, so the directory either does not exist, or we are not inside it
then
    #echo "Currently not in the directory current_profile"
    if [ -d $dirname ] 
    then
        printf "\nDirectory already exists.\nCalculation parameters:\n"
        grep "bond=\|fixed=\|in=\|out=\|up=\|down=\|MF=" $dirname/bond_integral.inp
        echo "Enter [y] to overwrite or any key to exit."; read accept
        if [ -z $accept ] || [ ! $accept == "y" ]  # if the variable is empty or different from "y", exit 
        then 
            exit
        else
            rm -f $dirname/*out $dirname/*inp # Remove leftover input and output to avoid confusion
        fi
    else
        mkdir $dirname
        cp MOL XDENS coord $dirname/
    fi  
fi

echo dista $atom1 $atom2
echo $(dista $atom1 $atom2)

distance=$( dista $atom1 $atom2 )
echo distance=$distance

# Coordinates of the centre of the bond:
BOND=$(centroid $atom1 $atom2)
Bx=$( echo $BOND | awk '{print $1} ' )
By=$( echo $BOND | awk '{print $2} ' )
Bz=$( echo $BOND | awk '{print $3} ' )

echo "Bond centre coordinates:"
printf "("$Bx"; "$By"; "$Bz")\n"

echo "Fixed point:"
fixed=1; read fixed # avoid empty value

in=2
printf "Define the starting point of integration\nPress [a] to enter atomic indices or press [v] the desired value\n"
read accept
echo $accept
if [ ! -z $accept ] && [ $accept == "a" ]
then
    echo "Centroid of atoms: "
    atoms=0; read atoms
    echo $atoms
	CENTROID=$(centroid $atoms)
	Dx=$( echo $CENTROID | awk '{print $1} ' )
	Dy=$( echo $CENTROID | awk '{print $2} ' )
	Dz=$( echo $CENTROID | awk '{print $3} ' )
	echo "Centroid:"
	printf "("$Dx"; "$Dy"; "$Dz")\n"
	awk -v Bx=$Bx -v By=$By -v Bz=$Bz -v Cx=$Cx -v Cy=$Cy -v Cz=$Cz BEGIN'{dist=sqrt( (Cx-Bx)*(Cx-Bx) + (Cy-By)*(Cy-By) + (Cz-Bz)*(Cz-Bz) ); print "Bond: ", Bx, By, Bz; print "Cent: ", Cx, Cy, Cz; print "diff: ", Cx-Bx, Cy-By, Cz-Bz; print dist}'
	in=$( awk -v Bx=$Bx -v By=$By -v Bz=$Bz -v Cx=$Cx -v Cy=$Cy -v Cz=$Cz BEGIN'{dist=sqrt( (Cx-Bx)*(Cx-Bx) + (Cy-By)*(Cy-By) + (Cz-Bz)*(Cz-Bz) ); print dist}' )
	
	printf  "Do you accept in="$in"?\nPress [n] to change\n"
	read accept;
	if [ ! -z $accept ] && [ $accept == "n" ]
	then
	    printf "in="; read in
	fi
    elif [ ! -z $accept ] && [ $accept == "v" ]
    then
	printf "in="; read in
	printf  "Do you accept in="$in"?\nPress [n] to change\n"
	read accept;
	if [ ! -z $accept ] && [ $accept == "n" ]
	then
	    printf "in="; read in
	fi
fi # end if accept in 

out=10
printf "Do you accept out=$out?\nPress [a] to enter atomic indices or press [v] the desired value\n"
read accept
echo $accept
if [ ! -z $accept ] && [ $accept == "a" ]
then
    echo "Centroid of atoms: "
    atoms=0; read atoms
	CENTROID=$(centroid $atoms)
	Dx=$( echo $CENTROID | awk '{print $1} ' )
	Dy=$( echo $CENTROID | awk '{print $2} ' )
	Dz=$( echo $CENTROID | awk '{print $3} ' )
	echo "Centroid:"
	printf "("$Dx"; "$Dy"; "$Dz")\n"
	out=$( awk -v Dx=$Dx -v Dy=$Dy -v Dz=$Dz -v Bx=$Bx -v By=$By -v Bz=$Bz BEGIN'{dist=sqrt( (Bx-Dx)*(Bx-Dx) + (By-Dy)*(By-Dy) + (Bz-Dz)*(Bz-Dz) ); print dist}' )
	
	printf  "Do you accept out="$out"?\nPress [n] to change\n"
	read accept;
	if [ ! -z $accept ] && [ $accept == "n" ]
	then
	    printf "out="; read out
	fi
    elif [ ! -z $accept ] && [ $accept == "v" ]
    then
	printf "out="; read out
	printf  "Do you accept out="$out"?\nPress [n] to change\n"
	read accept;
	if [ ! -z $accept ] && [ $accept == "n" ]
	then
	    printf "out="; read out
	fi
fi # end if accept out = 10


up=10
down=0
printf  "Do you accept up=10 and down=0?\nPress [n] to change\n"
read accept;
if [ ! -z $accept ] && [ $accept == "n" ]
then
    printf "up="; read up
    printf "down="; read down
fi

MF=z
printf  "Magnetic field direction = $MF. Do you accept?\nPress [n] to change\n"
read accept;
if [ ! -z $accept ] && [ $accept == "n" ]
then
    printf "MF="; read MF
fi

spacing=0.02
printf  "Do you accept spacing=[$spacing, $spacing, $spacing]?\nPress [n] to change\n"
read accept;
if [ ! -z $accept ] && [ $accept == "n" ]
then
    printf "spacing="; read spacing
fi


string="s/@bond@/$bond/; s/@fixed@/$fixed/; s/@distance@/$distance/; s/@up@/$up/; s/@down@/$down/; s/@MF@/$MF/; s/@spacing@/$spacing/g; s/@in@/$in/; s/@out@/$out/"
sed "$string" /home/mariavd/scripts/gimic/bond_integral.inp > ./$dirname/bond_integral.inp

printf "Performing dryrun...\n\n"
(cd ./$dirname/ && gimic --dryrun bond_integral.inp | grep "grid points" )
printf "\n\n"

echo "Ready for analysis"

