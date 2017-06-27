#!/usr/bin/env bash

function runtest_2D() {
    if [ $verbose -eq 1 ]
    then
	printf "\n\nPerforming test on $testname 2D current density\n\n"
    fi

    mkdir ../tmp/$testname
    cp ./$testname/MOL ./$testname/XDENS ../tmp
    cp ./$testname/2D/gimic.inp ../tmp/$testname
    (cd ../tmp/$testname && $gimicdir/gimic gimic.inp > gimic.test.out )


    # variable to track the number of the test executed
    i=0

    for file in "jvec.txt" "jvec.vti" "acid.txt" "jmod.txt"
    do
	i=$(( $i + 1 ))
	if diff ../tmp/$testname/$file ./$testname/2D/$file.ref >/dev/null
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

    rm -rf ../tmp/MOL rm -rf ../tmp/XDENS
} # end function runtest_2D


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

# Make a temporary directory for the test
mkdir ../tmp

#molecules="benzene C4H4"
molecules="benzene"

for testname in $molecules
do
    runtest_2D
done

if [ $verbose -eq 1 ]
then
    printf "\nSuccess of all tests:\n"
fi

rm -rf ../tmp

echo $success
