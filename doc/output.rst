

Description of generated files
==============================

* ``MOL``: contains molecular coordinates and basis functions, needed as input file for a gimic calculation
* ``XDENS``: contains information about the AO density and first order perturbed AO density, needed as input for a gimic calculation
* ``gimic.inp``: gimic input file, here either the cubic grid or a grid across a bond is specified, needed as input information for a gimic calculation
* ``gimic.out``: gimic output file, contains information about the integrated current strength susceptibility
* ``grid.xyz``: contains molecular coordinates in Ångstrøm and cube points (8) or integration plane corner points (4)
* ``mol.xyz``: contains molecular coordinates in Ångstrøm
* ``acid.txt``: contains information about a grid point and the ACID function value in atomic units at that point. x, y, z, f_acid(x,y,z), x,y,z are given in bohr
* ``acid.cube``: contains information about ACID function as Gaussian cube file
* ``acid.vti``: contains information about ACID function as VTK file compatible with ParaView
* ``jvec.txt``: contains information about a grid point x,y,z and the current density vector in that point jx, jy, jz (x,y,z, f_current(x,y,z), the point is given in bohr (it used to be in Ångstrøm)
* ``jvec.vti``: contains information about the current density vector function as VTK file compatible with ParaView
* ``jmod.txt``: contains a grid point and the modulus of the current density vector function, x,y,z, f_mod(x,y,z)
* ``jmod.vti``: contains the signed modulus of the current density vector function and can be visualized via two isosurfaces in ParaView
* ``jmod.cube``: contains the modulus of the current density vector function as Gaussian cube file
* ``jmod_quasi.cube``: contains the negative part of the signed modulus as Gaussian cube file
