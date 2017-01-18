#!/bin/bash

function runtest_3D() {
    if [ $verbose -eq 1 ]
    then
	printf "\n\nPerforming test on $testname 3D current density\n\n"
    fi

    (cd ./$testname/3D && $gimicdir/gimic gimic.inp > gimic.test.out )


    # variable to track the number of the test executed
    i=0

    for file in "jvec.txt" "jvec.vti" "acid.cube" "acid.txt" "acid.vti" "jmod.cube" "jmod_quasi.cube" "jmod.vti" "jmod.txt"
    do
	i=$(( $i + 1 ))
	if diff ./$testname/3D/$file ./$testname/3D/$file.ref >/dev/null
	then
#	    test$i=0
	    if [ $verbose -eq 1 ]
	    then
		echo test$i: $file Same
	    fi
	else
#	    test$i=1
	    success=$(( $success + 1 ))
	    if [ $verbose -eq 1 ]
	    then
		echo test$i: $file Different
	    fi
	fi
#	success=$(( $success + $(test$i) ))
    done

    if [ $verbose -eq 1 ]
    then
	echo Outcome of the test on $testname
	echo $success # successful result is success=0
    fi

}

arg="$2"
#echo Argument 2: $arg

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

#echo Verbose? $verbose

gimicdir="$1"

if [ -z $gimicdir ]
then
    echo "Gimic directory not specified as a command line argument"
    exit
fi

# Initialize the variable to check the success of the test runs
success=0

#molecules="benzene C4H4"
molecules="benzene"

for testname in $molecules
do
    runtest_3D
done

if [ $verbose -eq 1 ]
then
    printf "\nSuccess of all tests:\n"
fi

echo $success
