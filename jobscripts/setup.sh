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



# Prepare the base structure of the current profile script for a local machine:
file1=current-profile-header
file2=current-profile-local-submit
SCRIPT_OUT="current-profile-local.sh.in"
if [ -e $SCRIPT_OUT ]
then
    AGE1=$(( $(date -r $file1 +%s) - $(date -r $SCRIPT_OUT +%s) ))
    AGE2=$(( $(date -r $file2 +%s) - $(date -r $SCRIPT_OUT +%s) ))
else
    AGE1=1
    AGE2=1
fi
if [ "$AGE1" -gt 0 ] || [ "$AGE2" -gt 0 ]
then 
    cat current-profile-header > current-profile-local.sh.in
    cat current-profile-local-submit >> current-profile-local.sh.in
fi


# Prepare the base structure  of the current profile script for cluster
file1=current-profile-header
file2=current-profile-cluster-submit
SCRIPT_OUT="current-profile-cluster.sh.in"
if [ -e $SCRIPT_OUT ]
then
    AGE1=$(( $(date -r $file1 +%s) - $(date -r $SCRIPT_OUT +%s) ))
    AGE2=$(( $(date -r $file2 +%s) - $(date -r $SCRIPT_OUT +%s) ))
else
    AGE1=1
    AGE2=1
fi
if [ "$AGE1" -gt 0 ] || [ "$AGE2" -gt 0 ]
then 
    cat current-profile-header > current-profile-cluster.sh.in
    cat current-profile-cluster-submit >> current-profile-cluster.sh.in
fi






# Prepare the batch job scripts:

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


# Prepare all other scripts

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


