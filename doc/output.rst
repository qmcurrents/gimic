

Description of generated files
==============================

Note "file.txt" means text file format and "file.vti" stands for VTK file format, which is compatible with ParaView. The text files start with a given grid point (x,y,z) followed by either information about the modulus of the current density f(x,y,z) or the components of the current density vector Jx,Jy,Jz.  

* ``MOL``: contains molecular coordinates and basis functions, needed as input file for a gimic calculation
* ``XDENS``: contains information about the AO density and first order perturbed AO density, needed as input for a gimic calculation
* ``gimic.inp``: gimic input file, here either the cubic grid or a grid across a bond is specified, needed as input information for a gimic calculation
* ``gimic.out``: gimic output file, contains information about the integrated current strength susceptibility
* ``grid.xyz``: contains molecular coordinates in Ångstrøm and cube points (8) or integration plane corner points (4)
* ``mol.xyz``: contains molecular coordinates in Ångstrøm

closed-shell calculation:

* ``acid.vti``: contains information about ACID function 
* ``jvec.txt``: contains information about a grid point x,y,z and the current density vector in that point jx, jy, jz (x,y,z, f_current(x,y,z), the point is given in bohr (it used to be in Ångstrøm)
* ``jvec.vti``: contains information about the current density vector function 
* ``jmod.txt``: contains a grid point and the modulus of the current density vector function, x,y,z, f_mod(x,y,z)
* ``jmod.vti``: contains the signed modulus of the current density vector function 

only for open-shell calculations:

* ``jmodalpha.txt``: contains the signed modulus of the alpha current density vector function 
* ``jmodalpha.vti``: contains the signed modulus of the alpha current density vector function 
* ``jmodbeta.txt``: contains the signed modulus of the beta current density vector function 
* ``jmodbeta.vti``: contains the signed modulus of the beta current density vector function 
* ``jmodspindens.txt``: contains the signed modulus of the difference (alpha - beta) current density vector function 
* ``jmodspindens.vti``: contains the signed modulus of the difference (alpha - beta) current density vector function 

* ``jvecalpha.txt``: contains information about the alpha current density vector function 
* ``jvecalpha.vti``: contains information about the alpha current density vector function 
* ``jvecbeta.txt``: contains information about the beta current density vector function 
* ``jvecbeta.vti``: contains information about the beta current density vector function 
* ``jvecspindens.txt``: contains information about the difference (alpha - beta) current density vector function 
* ``jvecspindens.vti``: contains information about the difference (alpha - beta) current density vector function 

