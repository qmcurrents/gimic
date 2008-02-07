#!/bin/sh
##
##   Written by Jonas Juselius <jonas@iki.fi>
##   extended by Michael Harding <harding@uni-mainz.de> 
##

###-------- SGE parameters ----------
#$ -S /bin/bash
#$ -N test -o test.out -e qsub.out
#$ -V -cwd -j y -m beas -notify
#$ -M xxx@yyy.de
#$ -l vf=1700M
#$ -l arch=lx26-amd64
#$ -pe mpi 1
##-------- end SGE parameters ----------

###-------- PBS parameters ----------
#PBS -N @name@
#PBS -o run.out 
#PBS -e qsub.out
##PBS -q @queue@
#PBS -lnodes=@nodes@
#PBS -lwalltime=@wall@
#PBS -lpmem=@mem@MB
#PBS @qopts@
#PBS -m e
###-------- end PBS parameters ----------

###-------------LSF----------------------
#BSUB -n 2             # request number of processors
#BSUB -J bench
#BSUB -o bench.%J.log
#BSUB -e bench.%J.err
#BSUB -q normal      # normal queue
#BSUB -W 23:59            # wall time requested
#BSUB -R "span[ptile=1]"
#BSUB -x
###------------end LSF ------------------


PROG=@PROG@
QUESYS=$PROG/lib/quesys.sh

if [ -f $QUESYS ]; then
    .   $QUESYS
else
    echo "$QUESYS not found!"
    exit 1
fi

# if NOTIFY is set an e-mail notification is sent to that adress at the start
# and end of the job
NOTIFY="@notify@"

###------ JOB SPECIFIC ENVIRONMENT --###

#export OMP_NUM_THREADS=4
PATH=".:$PATH:$PROG/bin"
###------ JOB SPECIFIC DEFNITIONS ------###

#stop_on_crash=no
# possible jobs types are: mpich, lam, scali, mvapich or serial
paratype="@paratype@" 

# a job id is automatically added to the workdir
workdir="@workdir@/$USER" 
global_workdisk="@global@"
outdir="@outdir@"

###--- JOB SPECIFICATION ---###
input="gimic.inp XDENS mol"
initialize_job

# distribute input files to all nodes
distribute $input

# exenodes executes non MPI commands on every node. If the '-all' flag is
# given it will execute the command for every allocated CPU on every node.
#exenodes mycomand files foobar

# run() takes the following options:
#   -para run parallel program
#   -serial/-seq run non-parallel program on all nodes
#   -driver/-drv start a dirver program

@jobcmd@

# gather files from all nodes. 
# gather() accepts the follwing flags:
#  -tag      append the nodename to every file 
#  -maxsize  maximum size (kB) of files to copy back
#  -recurse  recurse through all subdirs
#  [list or pattern] at end to selectively copy files
#
#  if $gather_skiplist is set, the files in that list will be skipped
#  unconditionally
gather -maxsize 100000

finalize_job
###------ END  ------###
# vim:syntax=sh:filetype=sh
