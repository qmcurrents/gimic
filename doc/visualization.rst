
Vector plots
------------

We moved from the old "jvec.txt" output file to VTK file 
format. The information needed for vector plots is in the
file "jve.vti" and can be visualized using 
ParaView (https://www.paraview.org/). 

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

Note, animations are also possible using ParaView. 
 
(Signed) modulus density plots
------------------------------

This is possible with ParaView using the file "jmod.vti".

ACID plots
------------------------------

This is possible with ParaView using the file "acid.vti".

