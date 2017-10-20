

Installation
============

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

  $ make test


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
  $ ./setup --cmake-options="-DPYTHON_LIBRARY=/global/hds/software/cpu/eb3/Python/2.7.12-foss-2016b/lib/libpython2.7.so -DPYTHON_INCLUDE_DIR=/global/hds/software/cpu/eb3/Python/2.7.12-foss-2016b/include"
  $ cd build
  $ make
  $ make install
