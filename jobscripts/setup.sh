#!/bin/bash
SCRIPTS_IN=$(ls *in)
SCRIPTS_DIR=$(pwd)

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
    fi
done

echo; echo "USEFUL ALIASES:"

echo ' export GIMIC_HOME='${SCRIPTS_DIR/%jobscripts}

echo

echo ' alias gim="$GIMIC_HOME/jobscripts/gimic-run.sh" '
echo ' alias gcurrent="$GIMIC_HOME/jobscripts/current-profile.sh" '
echo ' alias gmol="python $GIMIC_HOME/jobscripts/turbo2gimic.py > MOL" '
echo ' alias grid="xmakemol -f grid.xyz &" '
echo ' alias dgrid="gimic --dryrun bond_integral.inp > /dev/null && xmakemol -f grid.xyz &" '
echo ' alias plotcurrent=" $GIMIC_HOME/jobscripts/plot-current-profile.sh" '
echo ' alias revplotcurrent="gnuplot $GIMIC_HOME/jobscripts/current-profile-plot-reverse.gnu" '
echo ' alias gpng="file=$(ls *para.eps) && convert -density 300 $file -resize 1024x1024 $file.png" '
echo ' alias gsq="$GIMIC_HOME/jobscripts/squares-profile.sh" '
echo ' alias geps="display *eps &" '
echo ' alias intprofile="$GIMIC_HOME/jobscripts/intprofile.sh" '
