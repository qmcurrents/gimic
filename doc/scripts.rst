

Interactive bash scripts
==========================

Introduction
----------------------------

GIMIC calculations can be prepared interactively using bash scripts. They are useful for adjusting the integration plane of choice, for starting current profile analysis, and submission of parallel tasks on a cluster with the SLURM schedulilng system (i.e. one where you submit jobs using ``sbatch``. There is a local version for the current profile analysis, which can run on a laptop, for example. Note that at the moment they are working with the ``coord`` file from TURBOMOLE. It will not be hard to make the input universal by using an XYZ file instead, so if you need it, let us know. The current density analyses employing GIMIC are most easily done if the magnetic field, pointing along the Z axis, is perpendicular to the molecular plane of a planar molecule. The orientation of the magnetic field can, of course, be modified accordingly. 

**NOTE:**  The script assumes that the ``MOL`` and ``XDENS`` files exist. They can prepared using ``turbo2gimic.py`` or ``Gaussian2gimic.py`` once the nuclear shielding calculation has finished. Also the coordinates need to be prepared in XYZ format, too. This can be done with the Turbomole script for conversion of the ``coord`` file to XYZ format.

::

   t2x coord > coord.xyz 


Setup 
-----------------------------
In the ``jobscripts`` directory there is a script called ``setup.sh``. It takes care of preparing the scripts for the specific machine. Remember, that it hardcodes the directories, so if you want to use the scripts on another computer, you need to execute ``./setup`` there. 

The scripts are still under development but the ones that are recommended for use at the moment are ``current-profile-cluster.sh``, ``current-profile-local.sh``, and ``gimic-run.sh``. The others are still in development phase, so if you are curious, you can try but functionality might be unclear or limited. 

There are some potentially useful bash aliases and functions in the file ``useful-aliases``, which you may want to add to your ``./bashrc``.

VMD is a highly-recommended tools for visualisation of molecules, in particular, the integration planes. 


Quick summary
-----------------------------

1. Choose an integration plane
2. If not obvious where to start integrating from, you can test its beginning and end points using the ``gimic-run.sh`` script and visualise the grid.xyz file. Adjust the gimic.inp file manually if necessary.
3. Run the current profile for the integration plane (``current-profile-local.sh`` if running on your laptop, or ``current-profile-cluster.sh`` if submitting as batch job on a cluster). 
4. Inspect the profile plot, check the sign of the currents, reverse the signs if necessary.
5. Make a better-quality plot using the ``plot-current-profile.sh`` script.


Detailed desctiption
-----------------------------

In the jobscripts directory you can find several bash scripts. The input for the scripts is based on selecting two atoms between which the integration plane lies. By default the integration plane is perpendicular to the bond between them. When you start the script, it asks for indices of the selected atoms according to the ``coord`` file. The atomic indices should be given in counter-clockwise manner. We have a sign problem, so for now we have to live with the fact that sometimes the diatropic current and the paratropic current can have the wrong sign but the magnitude is obtained correctly. Meaning that *it is not necessary to restart the calculation if the signs are inverted*. It can be post-processed. In planar molecules the choice of the two atoms is straight-forward but in non-planar ones, you need to choose the two atoms to be with very similar Z coordinates when placing the magnetic field in the Z direction. It is a known issue, which will hopefully be resolved in the future. 

**TIP:**   Check the indices of the atoms from the xyz file, for example with *xmakemol*. 

**NOTE:**   Atom indices start from 0 in VMD, whereas in GIMIC and xmakemol, they start from 1. 

A directory with a name include the atomic indices of the bond is created, for example ``current_profile_39.119``. This means that several current profiles can be analysed at the same time, each written to a directory of its own. If a directory with the same calculation already exists, it will ask for confirmation to overwrite it.

Next, you need to define the dimensions of the integration plane. It crosses the bond at the midpoint between the two attoms by default, however if needed, one can shift it to be closer to the one or the other atom. Usually we are interested in starting the integration from the centre of a ring to infinity, or maybe from the centre of the ring to another. Two parameters take care of that: the **start** value and the **out** value. These parameters can be defined either by manually entering a distance in atomic units, or using the atomic indices as defined in the ``coord`` file. By default, the script offers integration to infinity, which for practical purposes means 10 bohr. When entering these atomic indices, please type them on the same row and then press enter, for example,

::

   1 2 3 4 5 [ENTER]
   
**TIP:**   The integration can be defined to start or end at the centre of any number of atoms, including single atoms. 

The script will calculate the distance between the bond and the geometric centre of the selected atoms. It is recommended to add extra 0.5 bohr because the origin of the vortex lies in the geometrical centre of the ring only for symmetric molecules. Adding the additional distance from the bond allows the origin of the vortex to be distinguished easily on the current profile. Note that the **start** value can be negative, meaning that integration occurs on the other side of the midpoint of the bond. The height of the integration plane is usually defined to be 10 bohr above and below the molecular plane (the keywords for that are **up** and **down**). 

The current profile analysis is based on making thin slices of the defined integration plane, and executing a separate GIMIC calculation for each slice. The default width of the slice is 0.02 bohr but it can be modified if necessary for some reason. The number of grid points for the Gaussian quadrature used for the numerical integration is defined using the **spacing** keywords for X, Y and Z directions. The default values have been tested and are recommended. 

The magnetic field direction needs to be specified either as manually entered Cartesian coordinates, where the default is the Z direction (0; 0; -1), or using the program maximise-projection by Lukas Wirz. The program will be added to the GIMIC repository when it is finalized. More manipulation of the integration plane is done using the **fixed coordinate** and the **rotation** of the plane. The fixed coordinate is the third coordinate, which defines the integration plane. Its exact usage can be found in the source code at ``src/fgimic/grid.f90``. The script calculates an estimate, however, it is not always satisfactory. Alternatively, three Cartesian coordinates for the fixed point can be entered. This parameter is among the most difficult concepts in GIMIC, so one needs to get a feeling for it.  

After these parameters are specified, a dry run is performed to check if there are enough grid points for the Gaussian quadrature. If there are at least 9x9x1, then the result will be reliable. If that part succeeds, the input files for each of the slices of the integration plane are created. After that the script asks if a visualisation of the integration plane should be done. Selecting this option calls calculations of the grid at the first and last slices and writes a ``grid.check.xyz`` file. The current profile script can be put to background using ``CTRL+Z`` and the grid file opened, for example, in VMD. A useful representation of the plane can be done using the procedure below, which has to be placed in the ``~/.vmdrc`` file. It actually draws the integration planes for all the opened molecule in VMD.  

::

   proc intplanes {} {
   
       set loadedMolecules [molinfo list]
   
   	foreach molid $loadedMolecules {
   
   	    mol showrep $molid 0 off
   		mol modselect 0 $molid "all not element X Be"
   		set xel [atomselect $molid "element X"]
   		set coords [$xel get index]
   
   		for {set i 0} {$i < 4} {incr i} {
   		    lassign $coords i1 i2 i3 i4
   		}
   
   	    set c1 [atomselect $molid "index $i1"]
   		set c2 [atomselect $molid "index $i2"]
   		set c3 [atomselect $molid "index $i3"]
   		set c4 [atomselect $molid "index $i4"]
   
   		lassign [$c1 get {x y z}] pos1
   		lassign [$c2 get {x y z}] pos2
   		lassign [$c3 get {x y z}] pos3
   		lassign [$c4 get {x y z}] pos4
   
   		draw color red
   		set LINEWIDTH 6
   		draw line $pos1 $pos2 width $LINEWIDTH
   		draw line $pos3 $pos4 width $LINEWIDTH
   		draw line $pos1 $pos3 width $LINEWIDTH
   		draw line $pos2 $pos4 width $LINEWIDTH
   
   		set posHalf1 [ vecscale 0.5 [ vecadd $pos1 $pos2 ] ]
   		set posHalf2 [ vecscale 0.5 [ vecadd $pos3 $pos4 ] ]
   
   		set LINEWIDTH 3 
   		draw line $posHalf1 $posHalf2 width $LINEWIDTH
   	}
   
        mol representation CPK 0.600000 0.300000 50.000000 50.000000
   	mol color Element
   	mol material Opaque
   	mol addrep $molid
   	mol modselect 1 $molid "all not element X Be"
   }


After closing VMD, the current profile script should be brought back using ``fg`` and pressing enter again. On a cluster it will ask how many of the slices should be calculated in parallel, and what is the batch job limit. With that done, the ``sbatch`` command will be executed and one needs to wait for it to finish. When the job finishes, in the current profile directory there will be the ``current_profile.dat`` file. It lists four columns: the first one is the distance along the integration plane, and the next are net current, diatropic contribution and paratropic contribution respectively. These data are plotted as EPS files in the current profile directory. 

**NOTE:** One should make sure that the net current far from the molecule is diatropic (positive by convention). If not, the sign should be reversed using is wrong. The following alias can be used:

::

    alias revcurrent="mv current_profile.dat current_profile.dat.1 && awk '{printf \"%.6f\t%.6f\t%.6f\t%.6f\n\", \$1, -\$2, -\$4, -\$3}' current_profile.dat.1 > current_profile.dat"


The current strength of different peaks on the current profile can be obtained using the function below. It takes two numbers as arguments - distance along the integration plane. It sums the values of the current in the slices between these distances and returns the net, diatropic and paratropic currents. The distances can be obtained from the output in the file ``profile-points.out``. In the first column of the three sections there is distance along the integration plane, at which the diatropic, paratropic or net current vanish. One can estimate from the current profile plot which points are interesting and then take their actual values from the file. The data about the points can also be obtained by calling the script ``crit_pts.sh`` from the directory of the current profile. 

::

   function anprofile() { awk -v lower= -v upper= '{ if (( >= lower) && ( <= upper)) { total+=; dia+=; para+=; } } END { printf(nTotal current: %fnDiatropic: %fnParatropic: %fnn, total, dia, para); } ' current_profile.dat ;  };

   # used as:

   $ anprofile 0.22 1.5

   Net current: 0.066627
   Diatropic: 0.134426
   Paratropic: -0.067799


Finally, the current profile plot can them be done anew with ``plot-current-profile.sh``. It is still under development and it might not be very user-friendly at the moment, so feel free to ask questions about it. Please let us know if you have any further questions, bugs and ambiguous parts.




Some tips and advice
-----------------------------

**TIP:**  A short way to preview the integration plane without doing the calculations are dry runs:

::

   function dryrun() { gimic --dryrun "$@" > /dev/null ; xmakemol -f grid.xyz;  };


Choosing integration planes can be tricky. One way to get a better feeling for the current densities in a molecule is to start with the 3D calculation of the current density and exploring it in Paraview. You can use the new ``3D-run.sh`` script. It is rather basic at the moment, and the input file needs to be inspected and the grid checked before starting the calculation. Its aim is to create a decently large box around the molecule. For molecules larger than 100 atoms it is reasonable to use spacing of 1 bohr, otherwise the calculation takes too long. In small systems 0.5 bohr is a good choice; less than that might be an overkill, unless one needs close-up views of currents. The calculation can only run in serial. Once the calculation starts, GIMIC gives a good estimate of how long it will take. If it is unreasonable, size of the box or the spacing should be adjusted. The 3D calculation will give the ``jvec.vti`` file. It also prints the ``mol.xyz`` file. Paraview cannot handle XYZ files at the moment, so they need to be converted to CML format first. The bash function below can be used. It employs openbabel. In case openbabel is missing, the ``mol.xyz`` file can be saved as ``mol.cml`` in Avogadro and the first line of the function commented out. 

::

   function molecule() {
   babel -ixyz mol.xyz -ocml mol.cml
   awk '{ {FS="\""}; {OFS="\""};
        if ($1 ~ "<atom id") {
            if ($5 ~ "spinMultiplicity")
                { print $1, $2, $3, $4, $5, $6, $7, $8/0.526, $9, $10/0.526, $11, $12/0.526, $13 }
            else  { print $1, $2, $3, $4, $5, $6/0.526, $7, $8/0.526, $9, $10/0.526, $11 }
            }
        else print $0; }' mol.cml > mol-bohr.cml

   # It takes an XYZ file as an argument:

   $ molecule mol.xyz

The provided Paraview state file ``3D-LIC.pvsm``, it will ask about the location of the ``jvec.vti`` and ``mol-bohr.cml`` files. After they are loaded, it should present the line integral convolution (LIC) representation. In the *Slice* filter changing the z component of the origin permits exploring the current densities vertically. The number of arrows illustrating the current direction is adjusted from the *Glyph* filter. In the *Masking* group the selected Glyph Mode is *Every Nth Point*. Change the stride according to your preference. The length of the arrows is adjusted from the Scaling group, the *Scaling Factor* value. 

The 3D visual inspection help identifying which current vortices are interesting and where the integration plane would cross as few other vortices as possible. 


