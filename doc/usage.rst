

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
the ``$gimic`` keyword. When the calculation has finished Turbomole
writes out two files called ``CAODENS`` (AO density information) and
``XCAODEND`` (perturbed density information). Then run in the same
directory the ``turbo2gimic.py`` script (distributed with GIMIC) 
to produce the ``MOL`` and ``XDENS`` files.

::
    $ turbo2gimic.py > MOL

Note, open-shell calculations are not supported.


Running QChem and FERMION++
---------------------------

The scripts/input needed can be found in ``/tools/qchem``. 

Written by J. Kussmann, University of Munich, jkupc@cup.uni-muenchen.de

Convert Output of Q-Chemor FermiONs++ for GIMIC (TURBOMOLE format)

For a list of options, type 'qc2tm -h'...

qc2tm -h:

USAGE: qc2tm -t <qchem or fermions> -qcout <output-file> -scr
             <scratch-directory> -s2c (opt.) -openshell (opt.)

Note, open-shell calculations are not supported.

Running LSDalton
---------------- 

The scripts/input needed can be found in ``/tools/lsdalton2gimic``. 

Written by C. Kumar, University of Oslo, chandan.kumar@kjemi.uio.no

".GIMIC" needs to be added in LSDALTON.INP 

Run in the same directory "python lsdalton2gimic.py"

Note, open-shell calculations are not supported.

Running Gaussian
---------------- 

The scripts/input needed can be found in ``/tools/g092gimic``. 

This script has been provided by Vincent Liegeois University Namur.
vincent.liegeois@unamur.be


The program is made of two parts: 

1) Gaussian2gimic.py which is the main program

2) BasisSet.py which is a module file containing the functions to read the basis set and to do the transformations from Spherical to Cartesian.
   This file just need to be put in the same directory as Gaussian2gimic.py

Gaussian2gimic.py use optionparser to sets its different options.
Therefore, « Gaussian2gimic.py -h » will give you the full description.

There is only two options, only one of them is required.
  The command line to run the program on a formatted checkpoint file from gaussian is the following:

Gaussian2gimic.py -i file.fchk

This will produce two files: XDENS and MOL files

The extra option is -t followed by a XDENS coming from turbo2gimic.py

Gaussian2gimic.py -i file.fchk -t XDENS_turbo

This will, in addition to creating XDENS and MOL files, print in the standard output (on the terminal) a comparison with the values obtained from gaussian and from turbomole.
This is intended as a way to control that both numbers are similar.


For example, for a calculation on paranitroaniline with HF/6-311G(2df,2pd), the maximum error on the density matrix is 5.8E-5 while the maximum errors on the perturbed density matrices are 2.1e-2, 2.6e-2, 5.0e-2.

For the same molecule but with HF/cc-pVTZ, the maximum errors are: 4.8e-5, 2.5e-2, 5.4e-2, 7.8e-2.

ATTENTION, to have these agreements with turbomole, one need to specify "int=NoBasisTransform » in the Gaussian NMR calculation in order to prevent Gaussian from transforming the generalized contraction basis sets.


At last, the MOL file produced by Gaussian2gimic is slightly different from the one obtained by turbo2gimic.`

Indeed, turbo2gimic gives the basis set exactly as obtained on https://bse.pnl.gov/bse/portal with « optimized general contractions » checked.
BUT, the coefficients are not normalized with  « optimized general contractions » checked. but are normalized with « optimized general contractions » unchecked.

As an example, the first atomic orbital for cc-pVTZ basis set for the C with  « optimized general contractions » checked consist of a contraction of 8 GTOs while there are 10 GTOs in the contraction if « optimized general contractions » is unchecked.
The coefficients given in the website and in turbo2gimic are optimized for the contraction of the 10 GTOs not the 8.
Gaussian and therefore Gaussian2gimic gives the coefficients that are normalized for the contraction of the 8 GTOs.

Example input for benzene:

::

    %Chk=benzeneg09.chk
    %mem=2000mb

    #p B3LYP/Def2TZVP SCF=Tight NMR=GIAO Int=NoBasisTransform IOp(10/33=2) 

    Benzene Gaussian NMR example

    0 1
    C    1.2049777911    0.6956942520    0.0000000000
    C    1.2049777911   -0.6956942520    0.0000000000
    C    0.0000000000   -1.3913885041    0.0000000000
    C   -1.2049777911   -0.6956942520    0.0000000000
    C   -1.2049777911    0.6956942520    0.0000000000
    C    0.0000000000    1.3913885041    0.0000000000
    H    2.1430161769    1.2372709666    0.0000000000
    H    2.1430161769   -1.2372709666    0.0000000000
    H    0.0000000000   -2.4745419332    0.0000000000
    H   -2.1430161769   -1.2372709666    0.0000000000
    H   -2.1430161769    1.2372709666    0.0000000000
    H    0.0000000000    2.4745419332    0.0000000000
    
Running Gaussian creates a file "benzeneg09.chk" 
You need to convert this "*.chk" file to a formatted "*.fchk" file. 

::

$ formchk file.chk file.fchk  

Then you can proceed as described above and generate the MOL and XDENS
files with:

::

$ Gaussian2gimic.py --input=benzeneg09.fchk

Note, for open-shell cases you need to add "gfprint pop=regular iop(10/33=2)"
and use the Gaussian "log" file instead of the "fchk" file. 

Example input for benzene triplet dication 

::

    %LindaWorkers=cib26-2
    %NProcShared=20
    %Chk=benzeneg09.chk
    %mem=2000mb

    #p POP=FULL GFPrint nosymmetry B3LYP/DEF2TZVP SCF=Tight NMR IOp(10/33=2)

    Benzene Gaussian NMR example triplet dication

    2 3
    C    1.2049777911    0.6956942520    0.0000000000
    C    1.2049777911   -0.6956942520    0.0000000000
    C    0.0000000000   -1.3913885041    0.0000000000
    C   -1.2049777911   -0.6956942520    0.0000000000
    C   -1.2049777911    0.6956942520    0.0000000000
    C    0.0000000000    1.3913885041    0.0000000000
    H    2.1430161769    1.2372709666    0.0000000000
    H    2.1430161769   -1.2372709666    0.0000000000
    H    0.0000000000   -2.4745419332    0.0000000000
    H   -2.1430161769   -1.2372709666    0.0000000000
    H   -2.1430161769    1.2372709666    0.0000000000
    H    0.0000000000    2.4745419332    0.0000000000

::

$ mv file.out > file.log
$ Gaussian2gimic.py --input=file.log 

For the present example a current strength susceptibility of 8.4 nA/T
was calculated. 

Note the keyword "Int=NoBasisTransform" is only needed for reproducing
Turbomole based results. "NMR=GIAO" is not necessarily needed since 
using GIAO's is the default in G09.  
It seems there is a small difference between the keywords "nosymmetry"
and "Symmetry=None". The latter should only be used if "nosymmetry"
generates an error. 


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
