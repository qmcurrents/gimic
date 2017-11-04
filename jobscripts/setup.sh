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
