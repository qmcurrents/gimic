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
file3=functions-def
SCRIPT_OUT="current-profile-local.sh.in"

if [ $SCRIPT_OUT -ot $file1 ] || [ $SCRIPT_OUT -ot $file2 ] || [ $SCRIPT_OUT -ot $file3 ]  # [ FILE1 -ot FILE2 ]  -> True if FILE1 is older than FILE2, or is FILE2 exists and FILE1 does not.
then 
    cat current-profile-header > current-profile-local.sh.in
    cat current-profile-local-submit >> current-profile-local.sh.in
fi


# Prepare the base structure  of the current profile script for cluster

file1=current-profile-header
file2=current-profile-cluster-submit
file3=functions-def
SCRIPT_OUT="current-profile-cluster.sh.in"

if [ $SCRIPT_OUT -ot $file1 ] || [ $SCRIPT_OUT -ot $file2 ] || [ $SCRIPT_OUT -ot $file3 ]
then 
    cat current-profile-header > current-profile-cluster.sh.in
    cat current-profile-cluster-submit >> current-profile-cluster.sh.in
fi


# Prepare the batch job scripts:

file=jobscript.IN
SCRIPT_OUT=jobscript
sedstring="s:@SCRIPTS_DIR@:$SCRIPTS_DIR:"
if [ $SCRIPT_OUT -ot $file ] 
then 
    echo; echo "REMEMBER TO CHANGE THE BATCH SCRIPT jobscript-header TO SUIT YOUR CLUSTER BEFORE SETUP"
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
    if [ $SCRIPT_OUT -ot $file ] 
    then 
        sed "$sedstring" $file > $SCRIPT_OUT
        echo "Created script $SCRIPT_OUT."
	chmod +x $SCRIPT_OUT
    fi
done

echo; echo "Path to GIMIC:"
echo 'export GIMIC_HOME='${SCRIPTS_DIR/%jobscripts}
