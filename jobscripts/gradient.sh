wrkdir=$( pwd )

echo "Working directory:"; echo $wrkdir; echo

filenum=$(ls $wrkdir/*inp | wc -l)

cat /dev/null > $wrkdir/paratropic.dat #delete if it already exists
cat /dev/null > $wrkdir/diatropic.dat
cat /dev/null > $wrkdir/current.dat

out=$(grep out= $wrkdir/gimic.0.inp | grep -o -E '[0-9.]+')
start=$(grep in= $wrkdir/gimic.0.inp | grep -o -E '[0-9.]+')
tmp=$( awk -v start=$start BEGIN'{ print -start }' )
start=$(echo $tmp) # this is necessary because nowadays width starts with a negative value; in used to be positive
delta=$( awk -v out=$out -v start=$start 'BEGIN{ value=outr-start; delta=(value<0?-value:value); print delta }' )

echo $out, $start, $delta, $filenum

echo "Calculating the gradient..."

for (( i=0; i<$filenum; i++ ))
do
#    echo "i = $i"
    grep -A 2 "Induced current" $wrkdir/gimic.$i.out | awk -v wrkdir=$wrkdir '{ dia=sprintf("%s/diatropic.dat",wrkdir); para=sprintf("%s/paratropic.dat",wrkdir); if (NR == 2) printf("% f\n", $5) >> dia; else if (NR == 3) printf("% f\n", $5) >> para; }'
    grep "Induced current (nA/T)" $wrkdir/gimic.$i.out | awk -v i=$i -v start=$start -v delta=$delta -v wrkdir=$wrkdir '{ out=sprintf("%s/current.dat",wrkdir); printf("%5.2f\t% f\n", i*delta,$5) >> out; }'
done

paste $wrkdir/current.dat $wrkdir/diatropic.dat $wrkdir/paratropic.dat > $wrkdir/current_profile.dat
rm -f $wrkdir/paratropic.dat $wrkdir/diatropic.dat $wrkdir/current.dat

printf "\nData saved in current_profile.dat\n\n"

