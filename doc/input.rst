

The GIMIC input file
====================

The GIMIC input file is parsed by the ``getkw`` python parser, which
defines a grammar system based on sections and keywords in a recursive
manner. The input is in principle line oriented, but lines may be 
continued using a pipe symbol ``|`` at the end of a line. Furthermore, 
blanks and tabs are insignificant, with the exception of strings. 
Lines may be commented until end-of-line with a hash sign (#).

Sections are delimited by curly brackets ``section{...}``, and may have 
a keyword argument enclosed in parentheses ``section(argument){...}``.

Keywords come in two different types; simple keywords consisting of
integers, reals or strings (enclosed in double quote marks ``"string"``), 
and array keywords. Array keywords are enclosed in square brackets and the 
elements – integers, reals or strings – are delimited by a comma, e.g. 
``keyword = [0,1,2]``.

Keywords
--------

The top level section defines a few global parameters:

calc=cdens,integral
    This keyword determines what is to be calculated, and in what order.
    The possible options are: ’cdens’ – calculate current densities,
    ’integral’ – integrate the current flow through a cut-plane. Each of these options
    have their own respective sections to specify options and grids.

title
    Useless keyword, but since every program with a bit of self respect
    has a title, GIMIC also has one…

basis=MOL
    Name of the MOL file (eg. MOL or mol or whatever). A relative path can also be given.

density=XDENS
    Name of the density file (eg. XDENS). A relative path can also be given.

debug=1
    Set debug level. The higher the number, the more useless output one
    gets.

openshell=false
    Open-shell calculation
    
magnet=[Bx,By,Bz] 
    Define the direction of the external magnetic field by its vector 
    components

magnet\_axis=z
    Specify the magnetic field along a defined axis. Valid
    options are: i,j,k or x,y,z or X. “i,j,k” are the directions of the
    basis vectors defining the integration plane.
    “x,y,z” are the absolute fixed laboratory axis. 
    Note that ``magnet\_axis=X`` is used to specify the magnetic field 
    along the direction which is orthogonal to the molecular plane, but 
    parallel to the integration plane.

Grid section
-----------------

Grid types
~~~~~~~~~~~~
The currently implemented grids are

even
    The grid points are spaced uniformly 

base
    Uniformly distributed grid points in 2D or 3D. Used in current-density 
    calculations

bond
    Defining an integration plane passing through a particular chemical bond. 
    Used for the integration of the current strength. 

Custom grids can also be employed. 

Definition 
~~~~~~~~~~~~

A grid can be defined in the following manner:


Grid(type) {
  ...
}


Bond grid
~~~~~~~~~~~

The keywords specific to the **bond** type grid employed in the integration of the 
strength of the current density are:

bond=[a,b]
    Define the integration plane to cross the bond betweens atoms *a* and *b*. 
    The indices are chosen according to the MOL file (``coord`` file from Turbomole 
    calculations). Alternatively, two Cartesian coordinates can be employed via the 
    keywords ``coord1`` and ``coord2``.
    
coord1 = [x1,y1,z1]
    and
coord2 = [x2,y2,z2]
    Specify a pair of Cartesian coordinates instead of atomic indices. 
    The integration plane will cross the line between the two points. Alternatively, 
    use atomic indices using the keyword ``bond=[a,b]``.

fixpoint=c
    Specify an atomic index to use as the third point defining the integration plane. 
    Alternatively, specify an arbitrary point using the keyword ``fixcoord``.

fixcoord=[x, y, z]
    Specify an arbitrary point in 3D as the third point defining the integration plane.
    Alternatively, specify an atomic index using the keyword ``fixpoint``.
    
distance=r
    Define the distance between the two atoms or two Cartesian coordinates 
    between which the integration plane will cross. 
    
height=[-a, b]
    Specify the distance *a* between the bottom vertex of the integration plane 
    and the bond using the number *-a*. Specify the distance *b* between the top 
    vertex of the integration plane and the bond using the number *b*.
    
width=[-a,b]
    Specify the lengths *a* and *b* on both sides of the chemical bond.

type=gauss                  
    Use Gauss distribution of grid points for the Gauss quadrature
    
gauss_order=9               
    The order of the Gauss quadrature

Base grid
~~~~~~~~~~~~~~

The keywords for the **base** type grid employed in current-density calculations are:

origin=[x, y, z]
    This keyword specified the Cartesian cooridnate of the bottom left 
    corner of the plane or cube used during the current-density calculation. 

ivec=[x, y, z] 
    The direction of the vertical basis vector for the plane or cube.

jvec=[x, y, z] 
    The direction of the horizontal basis vector for the plane or cube. 
    The third vector specifying the cube is calculated as k = i × j, 
    therefore it is not given explicitly. 

length=[a, b, c] 
    The length of each side of the plane or cube. 

Universal keywords
~~~~~~~~~~~~~~~~~~~~

The following keywords are valid for both base and bond grids.

spacing=[a, b, c] 
    The distance between the grid points in the three directions (the basis 
    vectors i, j, and k). One should specify either ``spacing`` 
    or ``grid_points``.

grid_points=[a, b, c] 
    The specific number of grid points in each direction. One should specify either 
    ``spacing`` or ``grid_points``.

rotation=[a, b, c] 
    The angles of rotation in space.

rotation_origin=[x, y, z] 
    The point in space around which to rotate. If not specified, the rotation is 
    done at the middle of the bond.

Advanced section
-----------------

The keywords are given in the section:

Advanced {
   ...
}

The available keywords are:

lip_order=5
    Polynomial order of the Lagrange Interpolation Polynomials.
    If a calculation has been preformed on a even spaced grid, generate a
    grid suitable for Gaussian integration by doing Lagrange interpolation


spherical=off
    Use spherical cartesians (i.e. 5d/7f/10g…). This is usually handled
    automagically. Experts only.

diamag=on
    Turn on/off diamagnetic contributions (experts only)

paramag=on
    Turn on/off paramagnetic contributions (experts only)

GIAO=on
    Turn on/off gauge including atomic orbtitals (experts only)

screening=on
    Use screening to speed up calculations

screen\_thrs=1.d-8
    Screening threshold

Section: Essential
----------------------

The section contains some more specific keywords. 

ACID calculations
~~~~~~~~~~~~~~~~~~~~
ACID calculations can be defined by setting the calculation type to 
``calc=cdens`` and the grid to type base: ``Grid(base) {...}``

acid=on
    Turn on/off ACID calculation. It can only be done in current-density 
    calculations with the ``calc=cdens`` keyword and the respective grid. 

Modulus of the current density
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Calculate the mod(J) integral, this is useful to verify that the actual
integration grid is sensible in “tricky” molecules.

jmod=off
    Unless necessary, it should be turned off to save computational time.

Current-density calculations
------------------------------

Current density calculations are specified using the keyword ``calc = cdens``. 

ACID calculations can be performed in the current-density cal

acid=on
    Turn on/off ACID calculation

The produced files can be visualised in ParaView:

* ``acid.vti``: contains information about ACID function as VTK file compatible with ParaView
* ``jvec.vti``: contains information about the current density vector function as VTK file compatible with ParaView
* ``jmod.vti``: contains the signed modulus of the current density vector function and can be visualized via two isosurfaces in ParaView

Calculate the mod(J) integral, this is useful to verify that the actual
integration grid is sensible in “tricky” molecules.

Integration of the current strength
-------------------------------------
The calculation type has to be set to ``calc=integral`` and the grid to type bond. 
``Grid(bond) {...}``


