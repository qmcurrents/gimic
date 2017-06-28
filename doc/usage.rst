

Usage
=====

To run GIMIC three files are needed:

#. A file containing the effective one-particle density, and the
   magnetically perturbed densities in AO basis

#. A MOL file, with information on molecular geometry and basis sets

#. A GIMIC input file

The following sections explains how to obtain the density file and the
MOL file using either ACES2 or Turbomole.

Running ACES2
-------------

Using ACES2, the special driver script ’\ ``xgimic2.sh``\ ’ must be used
to run the NMR shielding calculation. Modify the script to suit your
needs (and set the paths correctly). If the NMR calculation is done with
symmetry, the MOL file must be converted to C1 symmetry using the script
``MOL2mol.sh``, prior to running GIMIC.

Example ZMAT:

::

    CO2
    O    2.14516685791074   0.00000000000000      0.00000000000000
    C    0.00000000000622   0.00000000000000      0.00000000000000
    O   -2.14516685791393   0.00000000000000      0.00000000000000

    *ACES2(CALC=CCSD,BASIS=tzp,UNITS=BOHR
    COORD=CARTESIAN
    MEMORY=250000000
    REFERENCE=RHF
    SYMMETRY=ON
    PROPERTY=NMR
    MULTIPLICITY=1
    CHARGE=0
    SCF_MAXCYC=200,CC_MAXCYC=150,CC_EXPORDER=40
    CC_CONV=10,SCF_CONV=10,LINEQ_CONV=10,CONV=10
    LINEQ_EXPAN=30)

Run ACES2 via xgimic2.sh to produce the XDENS file:

::

    $ xgimic2.sh --cc >aces2.out &

Convert the symmetry adapted MOL file to C1 symmetry:

::

    $ MOL2mol.sh

The new MOL file is now called mol.

Running Turbomole
-----------------

Starting with Turbomole 5.10, the GIMIC interface is part of the
official distribution. To produce the necessary files to run GIMIC, you
first need to optimize the wavefunction/density of the molecule, before
running the ``mpshift`` program to produce the perturbed densities.
Before you run ``mpshift`` you need to edit the ``control`` file and add
the ``$gimic`` keyword. When the calculation has finished run the
``turbo2gimic.py`` script (distributed with GIMIC) to produce the
``mol`` and ``XDENS`` files.

Running GIMIC
-------------

To run gimic you need to have at least three files: The gimic input file
(gimic.inp), the compound density file (XDENS) and the compound basis
set and structure file (mol). Copy the example gimic.inp (in the
``examples/`` directory) to your work directory, edit to your needs, and
execute

::

    $ gimic [--mpi] [gimic.inp] >gimic.out

To produce the mol and XDENS files:

- CFOUR: Do a normal NMR calculation and then run the 'xcpdens' program
  distributed with GIMIC to make the XDENS file. Then run the MOL2mol.sh
  script to produce the mol file.

- Turbomole: Add the $gimic keyword to the control file and then run mpshift
  as normal to produce a XDENS file. Then run the turbo2mol.py script to
  create the mol file from the coord and basis files.

Before doing the actual calculation it might be a good idea to check
that the grids are correct, run:

::

    $ gimic --dryrun

and examine the .xyz files that GIMIC produces. If they look ok, simply
run

::

    $ gimic

If you want to run the parallel version, there is a wrapper script
called ’\ ``qgimic``\ ’ (see ``qgimic –help`` for a list of command line
options) to produce a generic run script for most queueing systems. Eg.
to set up a parallel calculation with 8 CPUs, 1 h time and 200 MB memory
to be run in ``/work/slask``

::

    $ qgimic -n 8 -t 01:00 -m 200 /work/slask

This produces a ’gimic.run’ file. Edit this file and make sure it’s ok,
and then submit it to the queueing system:

::

    $ qsub gimic.run
