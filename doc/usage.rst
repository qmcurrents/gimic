******
Usage
******

Three files are needed in order to run  GIMIC:

#. The ``XDENS`` file containing the effective one-particle density, and the
   magnetically perturbed densities in AO basis;
   
#. The ``MOL`` file with information on molecular geometry and basis sets;

#. The GIMIC input file. The default name is ``gimic.inp``, however, any other 
   file name is acceptable.
   

The ``XDENS`` and ``MOL`` files are obtained externally using other quantum 
chemistry programs which support an interface to gimic. Details are 
presented below. 

Before doing the actual calculation with GIMIC, it might be a good 
idea to check if the grid is correctly defined. The ``dryrun`` flag
produces the ``grid.xyz`` file which can be examined in any molecular
viewer.

::

    $ gimic --dryrun

The calculation with GIMIC is started by:

::

    $ gimic [gimic.inp] > gimic.out

The square brackets mean that it is not necessary to give the name
of the input file if it is called ``gimic.inp`` since it is recognised 
by default. 

Interfaces to GIMIC
####################

The following sections explains how to obtain the ``XDENS`` and the ``MOL`` files.

Running TURBOMOLE
*******************

As of TURBOMOLE 5.10, the GIMIC interface is part of the
official distribution. To produce the necessary files to run GIMIC, you
first need to optimize the wavefunction/density of the molecule, before
running the ``mpshift`` program which produces the perturbed density matrices.
Before you run ``mpshift``, you need to edit the ``control`` file and add
the ``$gimic`` keyword. When the calculation has finished, TURBOMOLE
writes two files called ``CAODENS`` (AO density data) and
``XCAODENS`` (perturbed density data). In case either of them is 
missing, the ``mpshift`` calculation has not converged. It is necessary to add the 
``$csmaxiter N`` keyword to the control file, where ``N`` is larger than 
the default value of 30. Afterwards, run the ``turbo2gimic.py`` script (distributed 
with GIMIC) in the same directory to produce the ``MOL`` and ``XDENS`` files.

::

    $ turbo2gimic.py > MOL


Running LSDalton
******************

The scripts/input needed can be found in ``/tools/lsdalton2gimic``. 

Written by C. Kumar, University of Oslo, chandan.kumar@kjemi.uio.no.

``.GIMIC`` needs to be added in ``LSDALTON.INP``.

Run in the same directory:

:: 

   python lsdalton2gimic.py

Running GAUSSIAN
*****************

The nuclear shielding calculation on Gaussian needs to be performed by
including the keyword ``IOp(10/33=2)`` in order to print the perturbed density matrices
in the output file. Explicitly specifying ``NMR=GIAO`` is not necessary since
using GIAO's is the default in G09. The keyword ``Int=NoBasisTransform`` is
needed in order to prevent Gaussian from transforming the generalized
contraction basis sets. It ensures that the results will match with the ones
obtained in Turbomole. 

The general structure of the input file looks like this:

:: 

    %chk=file.chk
    #P <method>/<basis> nmr pop=regular Int=NoBasisTransform IOp(10/33=2)
    <empty line> 
    title 
    <empty line> 
    <charge> <multiplicity>
    <coordinates> 
    <empty line>


In the above, the lines with the ``<`` ``>`` symbols are supposed to be modified. These bracket symbols are not part of the actual input file. 

If the basis set needs to be specified explicitly, the input file is structured as follows:

::

    %chk=file.chk
    #P <method> nmr pop=regular Int=NoBasisTransform IOp(10/33=2)
    <empty line> 
    title 
    <empty line> 
    <charge> <multiplicity>
    <coordinates> 
    <empty line>
    <basis set specification>
    <****>
    <empty line> 

If ECPs are needed, then the method specification line should like as in the example below.

:: 

    #P TPSSTPSS/GenEcp nmr pop=regular Int=NoBasisTransform IOp(10/33=2)


It seems there is a small difference between the keywords ``nosymmetry``
and ``Symmetry=None``. The latter should only be used if ``nosymmetry``
generates an error. 

The scripts/input needed can be found in ``/tools/g092gimic``. 

This script has been provided by Vincent Liegeois from the University of Namur,
vincent.liegeois@unamur.be.

The tool consists of two parts: 

1) ``Gaussian2gimic.py`` which is the main script

2) ``BasisSet.py`` which is a module file containing the functions to read the basis set and to do the transformations from Spherical to Cartesian.
   This file just needs to be put in the same directory as ``Gaussian2gimic.py``.

``Gaussian2gimic.py`` uses optionparser to sets its different options.
Therefore, ``Gaussian2gimic.py -h`` will give you the full description.

There are two options: ``-i`` and ``-t``, however the latter is optional.
The command which runs the script on a formatted checkpoint file from GAUSSIAN is the following:

:: 
    
    Gaussian2gimic.py -i file.fchk

It will produce the ``XDENS`` and ``MOL`` files

The extra option ``-t`` accept the argument of the ``XDENS`` file, which matches the output of the interface to TURBOMOLE obtained with  ``turbo2gimic.py``.

:: 
    
    Gaussian2gimic.py -i file.fchk -t XDENS_turbo

In addition to creating ``XDENS`` and ``MOL`` files, the script prints to the terminal a comparison with the values obtained from GAUSSIAN and from TURBOMOLE.
This is intended as a way to make sure that both numbers are similar.

For example, for a calculation on paranitroaniline with HF/6-311G(2df,2pd), the maximum error on the density matrix is 5.8E-5 while the maximum errors on the perturbed density matrices are 2.1e-2, 2.6e-2, 5.0e-2.

For the same molecule but with HF/cc-pVTZ, the maximum errors are: 4.8e-5, 2.5e-2, 5.4e-2, 7.8e-2.

The ``MOL`` file produced by ``Gaussian2gimic`` is slightly different from the one obtained by turbo2gimic.`

Indeed, ``turbo2gimic.py`` gives the basis set exactly as obtained on https://bse.pnl.gov/bse/portal with "optimized general contractions" checked.
However, the coefficients are not normalized with  "optimized general contractions" checked. They are normalized when "optimized general contractions" is unchecked.

For example, the first atomic orbital in the cc-pVTZ basis set for carbon with  "optimized general contractions" checked consist of a contraction of 8 GTOs, while there are 10 GTOs in the contraction if "optimized general contractions" is unchecked.
The coefficients given in the website and in ``turbo2gimic.py`` are optimized for the contraction of the 10 GTOs rather than 8.
Gaussian and therefore ``Gaussian2gimic.py`` give the coefficients that are normalized for the contraction of the 8 GTOs.

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
    
Running Gaussian creates a file ``benzeneg09.chk``
You need to convert this ``*.chk`` file to a formatted ``*.fchk`` file. 

::

    $ formchk file.chk file.fchk  

Then you can proceed as described above and generate the ``MOL`` and ``XDENS``
files with:

::

    $ Gaussian2gimic.py --input=benzeneg09.fchk

Note that for open-shell cases you need to add "gfprint pop=regular iop(10/33=2)"
and use the Gaussian ``*.log`` file instead of the ``*.fchk`` file. 

Example input for the triplet dication of benzene:

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

When the calculation completes, run in the terminal:

::

   mv file.out > file.log
   Gaussian2gimic.py --input=file.log 

For the present example, a current strength susceptibility of 8.4 nA/T
was calculated. 


Running QChem and FERMION++
****************************

The scripts/input needed can be found in ``/tools/qchem``. 

Written by J. Kussmann, University of Munich, jkupc@cup.uni-muenchen.de.

Convert Output of Q-Chemor FermiONs++ for GIMIC (TURBOMOLE format)

For a list of options, type:

:: 

    qc2tm -h
    
which prints

::

    USAGE: qc2tm -t <qchem or fermions> -qcout <output-file> -scr
             <scratch-directory> -s2c (opt.) -openshell (opt.)

Running CFOUR
**************

Do a normal NMR calculation and then run the ``xcpdens`` program
distributed with GIMIC to make the ``XDENS`` file. Then run the ``MOL2mol.sh``
script to produce the ``MOL`` file.
  
Running ACES2
**************

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

Run ACES2 via ``xgimic2.sh`` to produce the ``XDENS`` file:

::

    $ xgimic2.sh --cc >aces2.out &

Convert the symmetry-adapted ``MOL`` file to C1 symmetry:

::

    $ MOL2mol.sh

The new ``MOL`` file is now called ``mol``.
