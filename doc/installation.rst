

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
