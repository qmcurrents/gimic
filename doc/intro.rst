

Introduction
============

This is the GIMIC program for calculating magnetically induced currents
in molecules. For this program produce any kind of useful information,
you need to provide it with an AO density matrix and three (effective)
magnetically perturbed AO density matrices in the proper format.
Currently only recent versions of ACES2 (CFOUR), Turbomole, QChem, 
LSDalton, FERMION++, Gaussian can produce these
matrices. Dalton is in the works.If you would like to add your
favourite program to the list please use the source, Luke.

-  For instructions how to compile and install this program refer to the
   INSTALL file in the top level directory.

-  There is an annotated example input in the examples/ directory.

-  For information on command line flags available run:
   ’\ ``gimic –help``\ ’

The following features have been implemented in the program

-  Current densities in 2D or 3D

-  The (signed) modulus of the current

-  The divergence of the current (this is useful for checking gauge
   invariance vs. gauge independence). Note, this module does not
work at the moment.

-  Vector representation of the current in 2D or 3D

-  Integration of the current flow through defined cut-planes in
   molecules

-  Open-shells and spin currents

-  Parallel execution through MPI (optional). Note, this works only if
you check out the "stable" branch.

Utility programs to extract the AO density and perturbed densities 
from ACES2, Turbomole, QChem, LSDALTON and FERMION++
calculations are included in the GIMIC source distribution, see tools
directory. Turbomole
5.10 and newer also has the GIMIC interface built in.
