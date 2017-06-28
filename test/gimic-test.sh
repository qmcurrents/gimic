#!/usr/bin/env bash

# radovan: this is to figure out the location of this script
# http://stackoverflow.com/questions/59895/can-a-bash-script-tell-what-directory-its-stored-in
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
    SCRIPT_DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
    SOURCE="$(readlink "$SOURCE")"
    # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
    [[ $SOURCE != /* ]] && SOURCE="$SCRIPT_DIR/$SOURCE"
done
SCRIPT_DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

cd $SCRIPT_DIR

function runtest() {
    if [ $verbose -eq 1 ]
    then
	printf "\n\nPerforming test on $testname bond integral\n\n"
    fi

    # Create a temporary directory for each test molecule, copy the input file there, then execute the calculation in it
    mkdir ../tmp/$testname

    cp ./$testname/int/gimic.inp ../tmp/$testname 
    cp ./$testname/MOL ./$testname/XDENS ../tmp
    (cd ../tmp/$testname  && $gimicdir/gimic gimic.inp > gimic.test.out )

    diatropic=$(grep -A 2 "Induced current" ../tmp/$testname/gimic.test.out | awk '{ if (NR == 2) printf("% f\n", $5); }')
    paratropic=$(grep -A 2 "Induced current" ../tmp/$testname/gimic.test.out | awk '{ if (NR == 3) printf("% f\n", $5); }')
    total=$(grep "Induced current (nA/T)" ../tmp/$testname/gimic.test.out | awk '{printf("% f\n", $5); }')

    # Debugging
    if [ $verbose -eq 1 ]
    then
	echo Calculated values:
	echo $diatropic $paratropic $total
    fi

    dia_ref=$(grep -A 2 "Induced current" ./$testname/int/gimic.out | awk '{ if (NR == 2) printf("% f\n", $5); }')
    para_ref=$(grep -A 2 "Induced current" ./$testname/int/gimic.out | awk '{ if (NR == 3) printf("% f\n", $5); }')
    total_ref=$(grep "Induced current (nA/T)" ./$testname/int/gimic.out | awk '{printf("% f\n", $5); }')

    # Debugging
    if [ $verbose -eq 1 ]
    then
	echo Reference values:
	echo $dia_ref $para_ref $total_ref
    fi

    # avoid returning success value by default
    test1=2
    test2=2
    test3=2

    test1=$( awk -v diatropic="$diatropic" -v dia_ref="$dia_ref" 'BEGIN{ if ((dia_ref - diatropic) < 1e-5) { print 0; } else {print 1; } }'  )
    test2=$( awk -v paratropic="$paratropic" -v para_ref="$para_ref" 'BEGIN{ if ((para_ref - paratropic) < 1e-5) { print 0; } else {print 1; } }' )
    test3=$( awk -v total="$total" -v total_ref="$total_ref" 'BEGIN{ if ((total_ref - total) < 1e-5) { print 0; } else {print 1; } }' )

    #Debugging
    if [ $verbose -eq 1 ]
    then
	echo "Tests for diatropic, paratropic and total current respectively:"
	echo $test1 $test2 $test3
    fi

    success=$(( $success + $test1 + $test2 + $test3 ))
    if [ $verbose -eq 1 ]
    then
	echo Outcome of the test on $testname
	echo $success # successful result is success=0
    fi

    rm -rf ../tmp/XDENS ../tmp/MOL
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

# Make a temporary directory for the test
mkdir ../tmp

molecules="benzene C4H4"

for testname in $molecules
do
    runtest
done

if [ $verbose -eq 1 ]
then
    printf "\nSuccess of all tests:\n"
fi

rm -rf ../tmp

echo $success
