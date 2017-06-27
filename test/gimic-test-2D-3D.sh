#!/usr/bin/env bash

function runtest() {
    dim=$1

    # Make a temporary directory for the test
    mkdir ../tmp

    if [ $verbose -eq 1 ]
    then
	printf "\n\nPerforming test on $testname $dim current density\n\n"
    fi

    mkdir ../tmp/$testname
    cp ./$testname/MOL ./$testname/XDENS ../tmp
    cp ./$testname/$dim/gimic.inp ../tmp/$testname
    (cd ../tmp/$testname && $gimicdir/gimic gimic.inp > gimic.test.out )


    # variable to track the number of the test executed
    i=0

    if [ $dim = "2D" ]; then
       files="jvec.txt jvec.vti acid.txt jmod.txt"
    fi
    if [ $dim = "3D" ]; then
       files="jvec.txt jvec.vti acid.txt jmod.txt acid.cube acid.vti jmod.cube jmod_quasi.cube jmod.vti"
    fi

    for file in $files
    do
	i=$(( $i + 1 ))
	if diff ../tmp/$testname/$file ./$testname/$dim/$file.ref >/dev/null
	then
	    if [ $verbose -eq 1 ]
	    then
		echo test$i: $file Same
	    fi
	else
	    success=$(( $success + 1 ))
	    if [ $verbose -eq 1 ]
	    then
		echo test$i: $file Different
	    fi
	fi
    done

    if [ $verbose -eq 1 ]
    then
	echo Outcome of the test on $testname
	echo $success # successful result is success=0
    fi

    rm -rf ../tmp
}


arg="$2"

if [ -z $arg ]
then
    arg=0 # verbose off
fi

if [ $arg = "-v" ]
then
    verbose=1 # verbose on
else
    verbose=0 # verbose off
fi

gimicdir="$1"

if [ -z $gimicdir ]
then
    echo "Gimic directory not specified as a command line argument"
    exit
fi

# Initialize the variable to check the success of the test runs
success=0

molecules="benzene"

for testname in $molecules
do
    runtest 2D
    runtest 3D
done

if [ $verbose -eq 1 ]
then
    printf "\nSuccess of all tests:\n"
fi

echo $success
