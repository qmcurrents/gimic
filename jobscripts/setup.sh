#!/bin/bash
SCRIPTS_IN=$(ls *in)
SCRIPTS_DIR=$(pwd)

checkMaxProj=$(command -v maximise_projection)
if [ ! -e $checkMaxProj ]
then
    echo "The program maximise_projection is not found."
    echo "Compile it and export its path to use the automatic orientation of the magnetic field."
#    exit
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

# Prepare the batch job scripts:
# if [ -e jobscript ]
# then
#     AGE=$(( $(date -r jobscript.IN +%s) - $(date -r jobscript +%s) ))
# else 
#     # jobscript does not exist, so create it
#     AGE=1
# fi
# echo $AGE
# if [ "$AGE" -gt 0 ] 
# then 
#     # jobscript.IN is newer compared to jobscript
#     sed "$sedstring" $file > $SCRIPT_OUT
#     echo "Created script $SCRIPT_OUT."
#     cat jobscript-header > jobscript
#     cat jobscript.tmp >> jobscript
#     rm -rf jobscript.tmp
#     echo; echo "REMEMBER TO CHANGE THE BATCH SCRIPT jobscript TO SUIT YOUR CLUSTER"
# fi


file=jobscript.IN
SCRIPT_OUT=$(echo ${file/.IN/} )
sedstring="s:@SCRIPTS_DIR@:$SCRIPTS_DIR:"
if [ -e $SCRIPT_OUT ]
then
    AGE=$(( $(date -r $file +%s) - $(date -r $SCRIPT_OUT +%s) ))
else
    AGE=1
fi
if [ "$AGE" -gt 0 ] 
then 
    echo; echo "REMEMBER TO CHANGE THE BATCH SCRIPT jobscript TO SUIT YOUR CLUSTER"
    cat jobscript-header > $SCRIPT_OUT
    sed "$sedstring" $file >> $SCRIPT_OUT
    echo "Created script $SCRIPT_OUT."
    chmod +x $SCRIPT_OUT
fi

