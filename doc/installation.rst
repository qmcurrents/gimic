

Installation
============

GIMIC is written in pure Fortran 90/95, and thus requires a good F95
compiler to compile. GIMIC has been compiled and tested to work with the
following compilers:

-  GNU gfortran (4.2)

-  g95 (0.9)

-  Intel ifort (please note that ifort 9.0 does not work if
   optimizations are enabled)

-  Portland pgf90

GIMIC uses the standard GNU autoconf generated ’configure’ scripts to
examine your system and pick sensible defaults for most build variables.
For a complete list of available configuration options run

::

    $ ./configure --help

It’s recommended that GIMIC is properly installed after compilation,
although not strictly necessary. Here /opt/gimic will be used as the
install path. To configure and build gimic without parallel capabilities
run

::

    $ ./configure --prefix=/opt/gimic
    $ make install

Then simply add /opt/gimic/bin to your path (or make the appropriate
links in /opt/bin), and you should be up and running.

If you are not happy with the default compiler picked up by configure
you can override the default by doing

::

    $ FC=myfavouritef90 ./configure --prefix=/opt/gimic

If you want to build the ACES2 interface (xcpdens) you first need to
build ACES2, and then run

::

    $ ./configure --prefix=/opt/gimic --with-aces2=/path/to/ACES2/lib

If configure can’t find the BLAS library (only needed for xcpdens), you
need to specify where to look for it add the following flag to
configure:

::

    --with-blas-dir=/path/to/lib

If this does not work for some reason you can specify exactly how to
link against BLAS on your system:

::

    --with-blas='-L/path/to/my/blas -lfooblas -lwhatever'

Good luck!
