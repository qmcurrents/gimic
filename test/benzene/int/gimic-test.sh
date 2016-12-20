#!/bin/bash

gimicdir="$1"

if [ -z $gimicdir ]
then
    echo "Gimic directory not specified"
    exit
fi

$gimicdir/gimic gimic.inp > gimic.test.out 

diatropic=$(grep -A 2 "Induced current" gimic.test.out | awk '{ if (NR == 2) printf("% f\n", $5); }')
paratropic=$(grep -A 2 "Induced current" gimic.test.out | awk '{ if (NR == 3) printf("% f\n", $5); }')
total=$(grep "Induced current (nA/T)" gimic.test.out | awk '{printf("% f\n", $5); }')

#echo $diatropic
#echo $paratropic
#echo $total

test1=$( awk -v diatropic="$diatropic" 'BEGIN{ if ((16.771417 - diatropic) < 1e-5) { print 0; } }'  )
test2=$( awk -v paratropic="$paratropic" 'BEGIN{ if ((-4.969283 - paratropic) < 1e-5) { print 0; } }' )
test3=$( awk -v total="$total" 'BEGIN{ if ((11.802134 - total) < 1e-5) { print 0; } }' )

if [ $test1 -eq 0 ] && [ $test2 -eq 0 ] && [ $test3 -eq 0 ]
then
    echo "success"
fi
