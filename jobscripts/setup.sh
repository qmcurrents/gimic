


SCRIPTS_IN=$(ls *in)
SCRIPTS_DIR=$(pwd)

for file in $SCRIPTS_IN
do
    SCRIPT_OUT=$(echo ${file/.in/} )
    sedstring="s:@SCRIPTS_DIR@:$SCRIPTS_DIR:"
    sed "$sedstring" $file > $SCRIPT_OUT
    echo "Created script $SCRIPT_OUT."
done
