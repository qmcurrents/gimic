

Installation
============

Before compiling GIMIC you need to make sure that you have installed
ideally all of the packages collected in
the ``requirements.txt`` file.
You need minimum the ones listed below. ::

* cython
* numpy
* runtest == 2.2.0

A convenient way to install the packages listed in ``requirements.txt``
is to install the 
Anaconda2 (https://www.anaconda.com/distribution/) Python distribution
first. 
Then you can simply install all the packages listed in the file
``requirements.txt`` by using:: 

  $ conda install name-of-package

GIMIC requires CMake to configure and build. CMake is invoked via a front-end script called ``setup``::

  $ ./setup
  $ cd build
  $ make
  $ make install

To see all available options, run::

  $ ./setup --help

Branch "master"::
GIMIC requires CMake to configure and build.::

  $ mkdir build
  $ cd build
  $ cmake ../
  $ make
  $ make install

Test the installation with::

  $ cd build
  $ make test

Note, some tests may require Valgrind and will fail if this
debugger is not available. However, this is no need to worry if all
other tests pass. 


Parallelization
---------------

OpenMP parallelization is available::

  $ ./setup --omp

MPI parallelization is in the works.


Installation on Stallo supercomputer
------------------------------------

::

  $ git clone git@github.com:qmcurrents/gimic.git
  $ cd gimic
  $ module load Python/2.7.12-foss-2016b
  $ module load CMake/3.7.2-foss-2016b
  $ virtualenv venv
  $ source venv/bin/activate
  $ pip install -r requirements.txt
  $ ./setup
  $ cd build
  $ make
  $ make install


Using BLAS1 and BLAS2 routines
------------------------------

With GNU compilers use::

  $ ./setup --blas

With Intel compilers and MKL use::

  $ ./setup --fc=ifort --cc=icc --cxx=icpc --cmake-options="-D ENABLE_MKL_FLAG=ON"


Installation on a Mac
------------------------------

The main problem with installing GIMIC on a Mac is setting the
correct paths and force the program to find all libraries. 
Here we share some recommendations but 
can not guaratee that GIMIC will for sure work on your Mac. GIMIC has
been developed on a Ubuntu/Linux operating system. 
Note, there are for sure more elegant ways to get GIMIC installed
on a Mac. If you figure them out, please share. 

*   Install git and cmake. 
*   Install Anaconda2 then you should be able to install most of the
    recommended packages listed in "requirements.txt". 
*   For installing the C++ (gcc) and Fortran (gfortran) compilers you
    can use Xcode, Brew or MacPorts. You just need to make sure that
    GIMIC is able to find them. A way to solve this is editing the
    ".bashrc" file. 
*   Check if you have the following paths in your ".bashrc". If not
    add them. 
        * export PATH=/Users/your_username:/Users/your_username/anaconda/bin:$PATH 
        * PATH="/Applications/CMake.app/Contents/bin":"$PATH"
        * export PATH=/Users/your_username/gimic/build/bin:$PATH
        * export PATH="/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin"







