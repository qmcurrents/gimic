.. role:: math(raw)
   :format: html latex


Grids
=====

There are two principal types of grids; the simple ’std (or base)’ grid,
which is defined by a pair of (orthogonal) basis vectors, and the ’bond’
grid which is mostly useful for defining cut-planes through bonds for
integration. There is also a third grid type ’file’, which specifies a
file containing gridpoints For the exact format of this file please
refer to the source in grid.f90. Furthermore there are two types of
grids, evenly spaced or with grid points distributed for Gauss-Legendre
or Gauss-Lobatto quadrature. This is specified with the
’type=even\|gauss\|lobatto’ keyword. When a quadrature grid is specified
the order of the quadrature must also be specified with the
’gauss\_order’ keyword. The number of grid points in each direction is
specified either explicitly using either of the array keywords
’grid\_points’ or ’spacing’. If the chosen grid is not a simple even
spaced grid, the actual number of grid points will be adjusted upwards
to fit the requirements of the chosen quadrature.

The shape of the grid can also be modified by the ’radius’ key, which
specifies a cutoff radius. This can be useful for integration. Sometimes
it’s practical to be able to specify a grid relative to a well know
starting point. The ’rotation’ keyword specifies Euler angles for
rotation according to the x->y->z convention. Note that the magnetic
field is not rotated, unless it is specified with ’magnet\_axis=i,j or
k’.

GIMIC automatically output a number of .xyz files containing dummy
points to show how the grids defined actually are laid out in space.

Basic grids
-----------

The ’std’ grid is defined by giving an ’origin’ and two orthogonal basis
vectors ’ivec’ and ’jvec’ which define a plane. The third axis is
determined from :math:`\vec k=\vec i\times\vec j`. The array ’lengths’
specifies the grid dimensions in each direction.

grid(std)
~~~~~~~~~

::

  type=even
  origin=-8.0, -8.0, 0.0
      Origin of grid
  
  ivec=1.0, 0.0, 0.0
      Basis vector i
  
  jvec=0.0, 1.0, 0.0
      Basis vector j ( k = i x j )
  
  lengths=16.0, 16.0, 0.0
      Lenthts of (i,j,k)
  
  spacing=0.5, 0.5, 0.5
      Spacing of points on grid (i,j,k)
  
  grid\_points=50, 50, 0
      Number of gridpoints on grid (i,j,k)
  
  rotation=0.0, 0.0, 0.0
      Rotation of (i,j,k) -> (i’,j’,k’) in degrees. Given as Euler angles
      in the x->y->z convention.

Bond grids
----------

The ’bond’ type grids define a plane through a bond, or any other
defined vector. The plane is orthogonal to the vector defining the bond.
The bond can be specified either by giving two atom indices,
’bond=[1,2]’, or by specifying a pair of coordinates, ’coord1’ and
’coord2’. The position of the grid between two atoms is determined by
the ’distance’ key, which specifies the distance from atom 1 towards
atom 2. For analysing dia- and paramagnetic contributions, the positive
direction of the bond is taken to be from atom 1 towards atom 2. Since
one vector is not enough to uniquely defining the coordinate system
(rotations around the bond are arbitrary), a fixpoint must be specified
using either the ’fixpoint’ atom index or the ’fixcoord’ keyword. This
triple of coordinates is also used to fix the direction of the magnetic
field when the ’magnet\_axis=T’ is used.

The shape and size of the bond grid is specified by the keywords 'width' and
'height'.  Each of the four components is relative to a point on the line
connecting the two reference atoms/coordinates.  Typically (but not
necessarily) the first component of both ranges is negative and the second
component is positive.

::
    width=[-1.5, 5.0]
    height=[-5.0, 5.0]

grid(std):
~~~~~~~~~~

type=gauss\|lobatto
    Use uneven distribution of grid points for quadrature

bond=1,2
    Atom indices for bond specification

fixpoint=5
    Atom index to use for fixing the magnetic field and grid orientation

coord1=0.0, 0.0, 3.14
    Coordinate of atom 1

coord2=0.0, 0.0, -3.14
    Coordinate of atom 2

fixcoord=0.0, 0.0, 0.0
    Fixation coordinate

distance=1.5
    Place grid ’distance’ between atoms 1 and atom 2

gauss\_order=9
    Order for Gauss quadrature

spacing=0.5, 0.5, 0.5
    Spacing of points on grid (i,j,k) (approximate)

grid\_points=50, 50, 0
    Number of grid points on grid (i,j,k) (approximate)

height=-4.0, 4.0
    Grid size relative to grid center

width=-1.0, 6.0
    Grid size relative to grid center

radius=3.0
    Create a round grid by cutting off at radius

rotation=0.0, 0.0, 0.0
    Rotation of (i,j,k) -> (i’,j’,k’) in degrees. Given as Euler angles
    in the x->y->z convention.

