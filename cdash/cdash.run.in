#!/bin/sh -l
#
# Example PBS run file for testing a branch with different 
# compilers and settings. Should be reimplemented in Python.
#
###-------- PBS parameters ----------
#PBS -N cdash.run
#PBS -o cdash.out 
#PBS -e qsub.out
#PBS -j oe
#PBS -lnodes=c24-15:ppn=8
#PBS -lwalltime=01:00:00
#PBS -lpmem=1000MB
###-------- end PBS parameters ----------

dashdir="@DASHBOARD_DIR@"

[ "x$site" = "x" ] && site=`hostname -s`
site="$site.$USER"

case $site in 
    stallo*)
    module load cmake
    module load boost
    module load eigen
    module load google-test
    module load valgrind
    export BOOST_ROOT=$BOOSTHOME
    export MPIEXEC_PREFLAGS="--mca btl self,tcp"
    toolchain_name=Intel
    ;;
esac

if [ "x$PBS_JOBID" != "x" ]; then
    $dashdir/run_dashboard.py --site=$site
    $dashdir/run_dashboard.py --site=$site --memcheck=on

    export CXX=mpic++
    export CC=mpicc
    export FC=mpif90
    $dashdir/run_dashboard.py --site=$site --build=Release --mpi=2 --omp=2

    $dashdir/run_dashboard.py --site=$site --toolchain_name=GNU \
        --cc=gcc --cxx=g++ --fc=gfortran --coverage --memcheck
else
    $dashdir/run_dashboard.py --site=$site
fi

# ------ END  ------
# vim:syntax=sh:filetype=sh
