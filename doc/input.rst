

The GIMIC input file
====================

| The GIMIC input file is parsed by the getkw input parser, which
  defines a grammar based on sections and keywords in a recursive
  manner. The input consists of sections containing keywords and/or
  other sections, and so on. The input is in principle line oriented,
  but lines may be continued using a ’
| ’ at the end of a line. Furthermore, blanks and tabs are
  insignificant, with the exception of strings. Lines may be commented
  until end-of-line with a hash sign (#).

Sections are delimited by an opening ’’ and closing ’’, and may have a
keyword argument enclosed between ’(’ and ’)’.

Keywords come in two different types; simple keywords consisting of
integers, reals or strings (enclosed in “ ”), and array keywords. Array
keywords are enclosed in ’[’ ’]’ and elements – integers, reals or
strings – are delimited by ’,’.

Keywords
--------

The top level section defines a few global parameters:

dryrun=off
    Don’t actually calculate anything. Good for tuning grids, etc. Can
    also be specified on the command line.

mpirun=off
    [boolean] Run in parallel mode.

title
    Useless keyword, but since every program with a bit of self respect
    has a title, GIMIC also has one…

basis=mol
    Name of the MOL file (eg. MOL or mol or whatever)

density=XDENS
    Name of the density file (eg. XDENS)

spherical=off
    Use spherical cartesians (i.e. 5d/7f/10g…). This is usually handled
    automagically. Experts only.

debug=1
    Set debug level. The higher the number, the more useless output one
    gets.

diamag=on
    Turn on/off diamagnetic contributions

paramag=on
    Turn on/off paramagnetic contributions

openshell=false
    Open-shell calculation

screening=off
    Use screening to speed up calculations

screen\_thrs=1.d-8
    Screening threshold

show\_up\_axis=true
    Mark the “up” axis in .xyz files

calc=cdens,…
    This keyword determines what is to be calculated, and in what order.
    Possible options are: ’cdens’ – calculate current densities,
    ’integrate’ – integrate the current flow through a cut-plane, ’divj’
    – calculate the divergence of the current. Each of these options
    have their own respective sections to specify options and grids.

The current density
-------------------

Section: cdens
~~~~~~~~~~~~~~

Name of output file containing the current tensors

Name of output file containing the current vectors

Vector which specifies the direction of the magnetic field.

magnet\_axis=z] Specify the magnetic field along a defined axis. Valid
options are: i,j,k or x,y,z or T. “i,j,k” are the directions of the
basis vectors defining the computational grid after any Euler rotation.
“x,y,z” are the absolute fixed laboratory axis. “T“ is used for
integration and specifies the direction which is orthogonal to the
molecular plane, but parallel to the integration plane.

Scaling factor for plotting purposes.

Annihilate the diamagnetic contribution to the current. Experts only.

Annihilate the paramagnetic contribution to the current. Experts only.

Grid to be used for calculating the currents. See the ”Grids“ section
for a description of how to specify grids.

Produce files suitable for plotting with ’gnuplot’ or ’gopenmol’

vector=JVEC
    File to contain the current vector field (gnuplot friendly)

modulus=JMOD
    File to contain the modulus of the current density (gnuplot
    friendly)

nvector=NJVEC
    File to contain the normalized current vector field (gnuplot
    friendly). Mostly useful for debugging purposes.

gopenmol=jmod.plt
    File to contain the current density in a gopenmol friendly format.

Integration
-----------

Section: integral
~~~~~~~~~~~~~~~~~

magnet\_axis=T] Specify the magnetic field along the direction which is
orthogonal to the molecular plane, but parallel to the integration
plane.

Vector which specifies the direction of the magnetic field.

Calculate the mod(J) integral, this is useful to verify that the actual
integration grid is sensible in “tricky” molecules.

Integrate the tensor components

If a calculation has been preformed on a even spaced grid, generate a
grid suitable for Gaussian integration by doing Lagrange interpolation

Polynomial order of the Lagrange Interpolation Polynomials

Grid to be used for calculating the currents. See the ”Grids“ section
for a description of how to specify grids.

The divergence of the current field
-----------------------------------

Subsection: divj
~~~~~~~~~~~~~~~~

Vector which specifies the direction of the magnetic field.

magnet\_axis=z] Specify the magnetic field along a defined axis. Valid
options are: i,j,k or x,y,z or T. “i,j,k” are the directions of the
basis vectors defining the computational grid after any Euler rotation.
“x,y,z” are the absolute fixed laboratory axis. “T“ is used for
integration and specifies the direction which is orthogonal to the
molecular plane, but parallel to the integration plane.

Filename of gOpenMol plot

Grid to be used for calculating the currents. See the ”Grids“ section
for a description of how to specify grids.

The electronic density
----------------------

The GIMIC program can also produce plots of the electronic density. This
code is very rudimentary currently, and cannot produce densities of
specific MOs or ranges of MOs.

Section: edens
~~~~~~~~~~~~~~

density=’EDENS’
    Filename which contains AO density. XDENS is fine usually.

density\_plot=’edens\_plt.txt’
    File name of density plot

gopemol=edens.plt
    Filename of gOpenMol plot

grid(std)
    [subsection] Grid to be used for calculating the currents. See the
    ”Grids“ section for a description of how to specify grids.
