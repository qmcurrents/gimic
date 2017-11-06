#!/bin/bash
SCRIPTS_IN=$(ls *in)
SCRIPTS_DIR=$(pwd)

checkMaxProj=$(command -v maximise_projection)
if [ ! -e $checkMaxProj ]
then
    echo "The program maximise_projection is not found."
    echo "Compile it and export its path."
    exit
fi

for file in $SCRIPTS_IN
do
    SCRIPT_OUT=$(echo ${file/.in/} )
    sedstring="s:@SCRIPTS_DIR@:$SCRIPTS_DIR:"
    if [ -e $SCRIPT_OUT ]
    then
        AGE=$(( $(date -r $file +%s) - $(date -r $SCRIPT_OUT +%s) ))
    else
        AGE=1
    fi
    if [ "$AGE" -gt 0 ] 
    then 
        sed "$sedstring" $file > $SCRIPT_OUT
        echo "Created script $SCRIPT_OUT."
	chmod +x $SCRIPT_OUT
    fi
done

echo; echo "USEFUL ALIASES AND FUNCTIONS:"; 
echo "_________________________________________"; echo

echo 'export GIMIC_HOME='${SCRIPTS_DIR/%jobscripts}

echo

echo 'python $GIMIC_HOME/install/bin/turbo2gimic.py > MOL; '
echo 'alias gim="$GIMIC_HOME/jobscripts/gimic-run.sh" '
echo 'alias gmol="python $GIMIC_HOME/jobscripts/turbo2gimic.py > MOL" '
echo 'alias grid="xmakemol -f grid.xyz &" '
echo 'alias dgrid="gimic --dryrun bond_integral.inp > /dev/null && xmakemol -f grid.xyz &" '
echo 'alias plotcurrent=" $GIMIC_HOME/jobscripts/plot-current-profile.sh" '
echo 'alias revplotcurrent="gnuplot $GIMIC_HOME/jobscripts/current-profile-plot-reverse.gnu" '
#echo 'alias gpng="file=$(ls *para.eps) && convert -density 300 $file -resize 1024x1024 $file.png" '
echo 'alias gsq="$GIMIC_HOME/jobscripts/squares-profile.sh" '
echo 'alias geps="display *eps &" '
echo 'alias intprofile="$GIMIC_HOME/jobscripts/intprofile.sh" '
echo
echo "Note: on a cluster replace the name of the script below to current-profile-cluster.sh"
echo 'alias gcurrent="$GIMIC_HOME/jobscripts/current-profile-local.sh" '
echo 
echo "function anprofile() { awk -v lower=$1 -v upper=$2 '{ if (($1 >= lower) && ($1 <= upper)) { total+=$2; dia+=$3; para+=$4; } } END { printf("\nTotal current: %f\nDiatropic: %f\nParatropic: %f\n\n", total, dia, para); } ' current_profile.dat ;  }; "
echo
echo 'function dryrun() { gimic --dryrun "$@" > /dev/null ; xmakemol -f grid.xyz;  }; '
