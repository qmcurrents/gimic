
Vector plots
------------

Vector plots can be done for example with gnuplot or Mayavi2. 
The information needed is on "jvec.txt". 

Vector plots are also feasible with ParaView, use "jvec.vti" for this.

Streamlines
-----------

Using the PyNgl library:

Under ``/tools/Visualization`` you find a python script that is
based on the streamline.py script by R. Bast, University of Tromsø,
Radovan.Bast@uit.no

You need to install PyNgl on your computer.
https://www.pyngl.ucar.edu/

Use the "plot" script in the same directory wher you have
the files "coord.xyz" and "jvec.txt". 

Using ParaView: (https://www.paraview.org/)

Follow the instructions under ``/examples/benzene/ParaView-README``

File needed: "jvec.vti"
Note, you need to generate your own molecule file, eg. in cml format,
which is readable in ParaView. You can use for example Avogadro
(https://avogadro.cc/) and read in mol.xyz and export it to mol.cml. 

Note, animations are also possible using ParaView. 
 
(Signed) modulus density plots
------------------------------

This is possible using the files "jmod.cube" and "jmod_quasi.cube" 
and for example Jmol (http://jmol.sourceforge.net/) or 
VMD (http://www.ks.uiuc.edu/Research/vmd/). 
Alternatively, this is also feasible 
with ParaView using "jmod.vti".

