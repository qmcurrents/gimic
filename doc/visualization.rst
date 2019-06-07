
Vector plots
------------

The information needed for vector plots is in the
file "jvec.{vti,vtu}" and can be visualised using 
ParaView (https://www.paraview.org/). 

For the visualisation of the current vector field we use two different file
formats.  'vti'/ImageData for grids that consist of orthogonal basis vectors in
two or three dimensions and evenly spaced grid points, and
'vtu'/UnstructuredGrid for arbitrary sets of points.  vti/vtu files consist of
nodes (grid points) and cells.  In vti files the node positions are defined by
a range and spacing in each dimension and the (cuboid) cells are defined
implicitly.  In vtu files, the grid point coordinates are given explicitly and
each cell is defined by a cell type (10 for tetrahedral) and the indices of the
defining nodes.   For each node we write down a vector (the induced current)
while the cells are not assigned a value for typical visualisations (ie, the
value can be set to 0.0).  However, cells need to be defined---paraview
interpolates the values given at grid points within each cell.

In case of vti files, all required information is given by the grid and no
further input is required.

In case of vtu files we require an external program to create a tetrahedral
mesh.  One choice for such a program is TetGen
(http://wias-berlin.de/software/index.jsp?id=TetGen&lang=1).  The workflow is:

1. obtain a grid of your choice.  The file format should be: no header entries,
each line contains xyz of a grid point.  That is, the number of lines equals
the number of grid points.

2. run 'grd2node.sh <grid file>' to obtain 'grid.node' which can be read by tetgen.

3. run 'tetgen grid.node' to obtain 'grid.1.ele'. This file contains one
quadruple of indices for each tetrahedral cell.

4. run your gimic job to obtain 'jvec.vtu'.  The generated 'grid.1.ele' is read
automatically.

Check also the section on the interpretation of GIMIC
results and the Youtube examples. 

Streamlines
-----------

Follow the instructions under ``/examples/benzene/ParaView-README``

File needed: "jvec.vti"
Note, you need to generate your own molecule file, eg. in cml format,
which is readable in ParaView. You can use for example Avogadro
(https://avogadro.cc/) and read in mol.xyz and export it to mol.cml. 
Then you need to convert the mol.cml file to bohr. 

Note that animations are also possible using ParaView. 
 
(Signed) modulus density plots
------------------------------

This is possible with ParaView using the file "jmod.vti".

ACID plots
------------------------------

This is possible with ParaView using the file "acid.vti".

