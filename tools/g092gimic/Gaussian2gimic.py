#! /usr/bin/env python
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import BasisSet
import os.path
from optparse import OptionParser
import numpy
import math
import re

Table=['Xx', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn']

class FCHK (object) :

    def __init__ (self, fchkname='gaussian.fchk'):
        self.filename = fchkname

    def get_property(self, property):
        """
        Read a property name "property" in the fchk file
        return a numpy vector
        """
        f = open(self.filename, 'r')
        lines = f.readlines()
        f.close()
        # find the block "property"
        first = 0
        last = 0
        cols = {"R":5, "I":6}
        numpytype = {"R":"d", "I":"int"}
        stype = "R"
        size = 0
        for i, l in enumerate(lines):
            if property in l:
                l = l[len(property):]
                if "N=" in l:
                    # read array of data
                    pattern = r"\d+"
                    size = int(re.findall(pattern, l)[0])
                    typepattern = r"   [IR]   "
                    try:
                        stype = re.findall(typepattern, l)[0].strip()
                    except IndexError:
                        raise Exception("The type of data in fchk is not recognized")
                    first = i+1
                    nlines = (size-1)//cols[stype]+1
                    last = first+nlines
                    break
                else:
                    # single data
                    typepattern = r"   [IR]   "
                    try:
                        stype = re.findall(typepattern, l)[0].strip()
                    except IndexError:
                        raise Exception("The type of data in fchk is not recognized")
                    if stype == 'I':
                        pattern = r"\d+"
                        return int(re.findall(pattern, l)[0])
                    elif stype == 'R':
                        pattern = r"-?\d+\.\d*[ED]?[+\-]?\d*"
                        return float(re.findall(pattern, l)[0])
        lines = lines[first:last]
        if len(lines) == 0:
            return None
        # read the data
        data = []
        for line in lines:
            data.extend([float (fl) for fl in line.split()])
        data = numpy.array(data, numpytype[stype])
        if data.size != size:
            raise Exception("The number of data recovered [%i] is not the same as the size written in the file [%i]"%(data.size, size))
        return data

    def get_density(self):
        """
        Get the non-perturbed density
        """
        data = self.get_property("Total SCF Density")
        density = None
        if isinstance(data, numpy.ndarray):
            nbasis = self.get_property("Number of basis functions")
            index = 0
            density = numpy.zeros((nbasis, nbasis))
            for i in range(nbasis):
                for j in range(i+1):
                    density[i, j] = data[index]
                    index += 1

            # fill the other half of the density matrix
            density = density+density.transpose()-numpy.diag(density.diagonal())
        return density

    def get_spindensity(self):
        """
        Get the spin density
        """
        data = self.get_property("Spin SCF Density")
        density = None
        if isinstance(data, numpy.ndarray):
            nbasis = self.get_property("Number of basis functions")
            index = 0
            density = numpy.zeros((nbasis, nbasis))
            for i in range(nbasis):
                for j in range(i+1):
                    density[i, j] = data[index]
                    index += 1

            # fill the other half of the density matrix
            density = density+density.transpose()-numpy.diag(density.diagonal())
        return density

    def get_density_magnetic(self):
        """
        Get the perturbed density by magnetic field
        """
        data = self.get_property("Magnetic Field P1 (GIAO)")
        density = None
        if isinstance(data, numpy.ndarray):
            nbasis =  self.get_property("Number of basis functions")
            index = 0
            density = numpy.zeros((3, nbasis, nbasis))
            for idir in range(3):
                for i in range(nbasis):
                    for j in range(i+1):
                        density[idir, i, j] = data[index]
                        index += 1

                # fill the other half of the density matrix
                density[idir] = density[idir]-density[idir].transpose()
        return density

    def get_basisset(self):
        """
        read the basis set
        return a list of AO objects
        """
        nshell = self.get_property("Number of contracted shells")
        shelltype = self.get_property("Shell types")
        primitive_pershell = self.get_property("Number of primitives per shell")
        shell2atom = self.get_property("Shell to atom map")
        primitive_exp = self.get_property("Primitive exponents")
        primitive_coeff = self.get_property("Contraction coefficients")
        primitive_coeff_p = self.get_property("P(S=P) Contraction coefficients")
        atomicnumbers = self.get_property("Atomic numbers")
        # basis set list
        basisset = []
        start = 0
        for ishell in range(nshell):
            nprim = primitive_pershell[ishell]
            # get the exponents and coeffs for all the primitive of the shell
            exponents = primitive_exp[start:start+nprim]
            coefficients = primitive_coeff[start:start+nprim]
            coefficients_p = None
            if primitive_coeff_p is not None:
                if numpy.nonzero(primitive_coeff_p[start:start+nprim])[0].size !=0 :
                    coefficients_p = primitive_coeff_p[start:start+nprim]
            iatom = shell2atom[ishell]
            an = int(atomicnumbers[iatom-1])
            type = int(shelltype[ishell])
            shell = BasisSet.SHELL(iatom, an, type, exponents=exponents, coefficients=coefficients, coefficients_p=coefficients_p)
            basisset.append(shell)
            start = start + nprim
        if len(basisset) == 0:
            basisset =  None
        else:
            basisset =  BasisSet.BasisSet(basisset)
#            basisset.check_contractedcoeff()
        return basisset

class Gaussian (object) :

    def __init__ (self, logname='gaussian.log',w=None) :
        self.filename = logname

    def get_density(self):
        """
        Get the density matrix for alpha and beta electrons
        Need pop=regular or full keyword
        """
        f = open(self.filename, 'r')
        lines = f.readlines()
        f.close()
        pattern = r"\w* Density Matrix"
        size = 0
        densities = []
        for iline, line in enumerate(lines):
            if "basis functions," in line:
                size = int(re.findall(r"\d+", line)[0])
            elif size > 0 and re.search(pattern, line) is not None:
                ntab = (size-1) // 5 + 1
                nlines = (size+1) * ntab - 5*ntab*(ntab-1)//2
                densities.append(self.read_integrals(string=lines[iline+1:iline+1+nlines], size=size, symmetric=True))
        if len(densities) == 0:
            return None
        densities = numpy.array(densities)
        if densities.shape[0] == 1:
            densities = densities[0]
        return densities

    def get_density_magnetic(self):
        """
        Get the perturbed density by magnetic field for alpha and beta electrons
        Need IOP(10/33=2) keyword
        """
        f = open(self.filename, 'r')
        lines = f.readlines()
        f.close()
        pattern = r"P1 \w+ \(asymm, AO basis\)"
        size = 0
        densities = []
        for iline, line in enumerate(lines):
            if "basis functions," in line:
                size=int(re.findall(r"\d+", line)[0])
            elif size > 0 and re.search(pattern, line) is not None:
                ntab = (size-1) // 5 + 1
                nlines = (size+1) * ntab - 5*ntab*(ntab-1)//2
                densities.append(self.read_integrals(string=lines[iline+1:iline+1+nlines], size=size, symmetric=False))
        if len(densities) == 0:
            return None
        return numpy.array(densities)


    def read_integrals(self, string, size, symmetric=None):
        """
        Get integrals of the form:
                1             2             3             4             5
      1  0.100000D+01
      2  0.219059D+00  0.100000D+01
      3  0.000000D+00  0.000000D+00  0.100000D+01
      4  0.000000D+00  0.000000D+00  0.000000D+00  0.100000D+01
      5  0.000000D+00  0.000000D+00  0.000000D+00  0.000000D+00  0.100000D+01
      6  0.184261D+00  0.812273D+00  0.000000D+00  0.000000D+00  0.000000D+00
        """
        # read the data (triangular matrix divided into different block)
        integral = numpy.zeros((size,size))
        iblock = 0
        while len(string) > 0:
            end = size+1-iblock*5 #skip the first line that are the label of the  columns
            block = string[1:end]
            string = string[end:]
            for irow, blockline in enumerate(block):
                pattern = r"-?\d+\.\d*[ED]?[+\-]?\d*"
                alist = re.findall(pattern, blockline)
                line=[float(fl.replace('D','E')) for fl in alist]
                integral[irow+iblock*5, iblock*5:iblock*5+len(line)] = line
            iblock = iblock+1

        if symmetric is None:
            return integral
        elif symmetric:
            # fill the other half of the matrix if symmetric = True
            integral = integral+integral.transpose()-numpy.diag(integral.diagonal())
        else:
            # fill the other half of the matrix if antisymmetric (symmetric=False)
            integral = integral-integral.transpose()

        return integral

    def get_basisset(self):
        """
        read the basis set given by gfprint in log file
        return a dictionary with the basis set information of each atom
        The key is the label
        The value is a list made of tuple for each shell
        Shell[0] is the type of orbital
        Shell[1] is a numpy array with coeff and exp for each primitive
                 One line per primitive with an exponent and the coefficients
        """
        f = open(self.filename, 'r')
        lines = f.readlines()
        f.close()
        start = 0
        end = 0
        # get the block with the basis set informations
        for i, l in enumerate(lines):
            if "AO basis set" in l:
                start=i+1
            elif (start != 0) and ((l.strip().split()[0] != "Atom") and (l.strip().split()[0][:2] != "0.")) : # End of the basis set
                end = i
                break
        lines = lines[start:end]
        if len(lines) == 0:
            return None
        basisset = []
        while len(lines) > 0 :
            # line should look like:
            # Atom C1       Shell     1 S   6     bf    1 -     1          1.000000000000          2.000000000000          8.000000000000
            info = lines[0].split()
            label = info[1]
            # get iatom and an from label
            iatom = int(re.findall(r"\d+",label)[0])
            symbol = re.findall(r"[a-zA-Z]+",label)[0]
            an = Table.index(symbol)
            nprim = int(info[5])
            mult = int(info[9]) - int(info[7])+1
            typeoforbital = BasisSet.mult2type[mult]
            coord = [float(x) for x in re.findall(r"[+-]?\d+\.\d*", lines[0])]
            coord = numpy.array(coord)
            # Define a block that correspond to a shell
            block = lines[1:1+nprim]
            lines = lines[1+nprim:]
            # read the block of coeff and exp
            data = []
            for blockline in block:
                aline = [float(fl.replace('D', 'E')) for fl in blockline.split()]
                data.append(aline)
            data = numpy.array(data)
            exponents = data[:, 0]
            coefficients = data[:, 1]
            coefficients_p = None
            if data.shape[-1] == 3:
                coefficients_p = data[:, 2]
            shell = BasisSet.SHELL(iatom, an, typeoforbital, exponents=exponents, coefficients=coefficients, coefficients_p=coefficients_p, coord=coord)
            basisset.append(shell)

        if len(basisset) == 0:
            basisset =  None
        else:
            basisset = BasisSet.BasisSet(basisset)
        return basisset

usage = "usage: %prog [options] "
parser = OptionParser(usage=usage)

parser.add_option("-i", "--inputfile", dest="inputfile",
                  help="Input log or fchk file from an NMR calculation. For log file, you need to have specified the following keywords: gfprint pop=regular iop(10/33=2). gfprint: to print the database. pop=regular or pop=full: to print the density matrix. iop(10/33=2) to print the perturbed density matrices. NOTE: for OPENSHELL systems, only the log file can be used (with the three keywords above).",
                  action="store", type="string",default="none")

parser.add_option("-t", "--turbomole", dest="turbomole",
                  help="Arrange the XDENS like with Turbomole instead of like with CFOUR. Default: False (like CFOUR)",
                  action="store_true", default=False)

parser.add_option("", "--XDENS", dest="XDENS",
                  help="XDENS file to compare with the generated one. Default: none",
                  action="store", type="string", default="none")

(options, args) = parser.parse_args()

#check that the user introduce  one input filename
if options.inputfile == "none":
    parser.error("No inputfile given")

if not os.path.isfile(options.inputfile):
    parser.error("The file ['%s'] does not exist" % options.inputfile)



# open the fchk or log file
if os.path.splitext(options.inputfile)[1] == ".fchk":
    file=FCHK(options.inputfile)
    # Read the density matrices
    # Only for closeshell. For opensell, use log file
    density = file.get_density()
    if density is None:
        raise Exception("Error while reading the density matrix.")
    spindensity = file.get_spindensity()
    if spindensity is not None:
        parser.error("Read the log file for openshell systems.\nUse pop=regular gfprint iop(10/33=2) keywords.")
    density_m = file.get_density_magnetic()
    if density_m is None:
        raise Exception("Error while reading the perturbed density matrices.\nCheck that NMR keyword has been used.")
    # create the 'densities' numpy array to collect all of them
    densities = numpy.zeros((1, 4) +density.shape)
    densities[0, 0] = density
    densities[0, 1:] = density_m
elif os.path.splitext(options.inputfile)[1] == ".log":
    file = Gaussian(options.inputfile)
    density = file.get_density()
    if density is None:
        raise Exception("Error while reading the density matrix.\nCheck that pop=regular or pop=full keyword has been used.")
    density_m = file.get_density_magnetic()
    if density_m is None:
        raise Exception("Error while reading the perturbed density matrices.\nCheck that iop(10/33=2) keyword has been used.")
    # create the 'densities' numpy array to collect all of them
    if density.ndim == 2:# close shell
        densities = numpy.zeros((1, 4) +density.shape)
        densities[0, 0] = density
        densities[0, 1:] = density_m
    else:#open shell
        k = density.shape[-1]
        densities = numpy.zeros((2, 4, k, k))
        densities[0, 0] = density[0]
        densities[1, 0] = density[1]
        density_m = density_m.reshape((2, 3, k, k))
        densities[0, 1:] = density_m[0]
        densities[1, 1:] = density_m[1]
else:
    parser.error("The inputfile is not a fchk or log file")


#read the basis set
basisset = file.get_basisset()
if basisset is None:
    parser.error("No basis set information found.\nIf you read a log file, check that gfprint keyword has been used.")

basisset = basisset.split_SP()

# prepare a transformation matrix from Spherical-harmonic to Cartesian orbitals
# create the block diagonal matrix of transformation
if options.turbomole:
    CartOrdering = "turbomole"
else:
    CartOrdering = "cfour"

if basisset.spherical:
    tmat = basisset.Cart2Spher(square=False, real=True, order=CartOrdering, normalized=False)
else:
    tmat = basisset.Cart2Cart(order1=CartOrdering, order2="gaussian", normalized1=False, normalized2=True)
if tmat is None:
    raise Exception("Error when creating the transformation matrix")

nbasisset = basisset.to_Cartesian()

# arrange the order of the basis function like CartOrdering
# S for all atoms, P for all atoms, D for all atoms, ...
aos = nbasisset.to_AOs(CartOrdering)
if options.turbomole:
    index = aos.index_sort_by_type()
else:
    index = aos.index_sort_by_atom()
aos = aos.basis_from_index(index)

if densities.shape[2] != tmat.shape[1]:
    raise Exception("The shape of the density matrix [%i] is not the same as the shape of the transformation matrix [%i]"%(densities.shape[2], tmat.shape[1]))
if len(index) != tmat.shape[0]:
    raise Exception("The length of the list of index [%i] is not the same as the shape of the transformation matrix [%i]"%(len(index), tmat.shape[0]))
tmat = tmat[index, :]#new, old
if options.turbomole:
    nbasisset = nbasisset.basis_from_index(nbasisset.index_sort_by_type())
else:
    nbasisset = nbasisset.basis_from_index(nbasisset.index_sort_by_atom())

# apply the transformation to 'densities' array
# the transformation reorganized the basis functions
# with the shells organized by types
# with the angular momenta organized like CartOrdering
nspin = densities.shape[0]
ndensities = numpy.zeros((nspin, 4, tmat.shape[0], tmat.shape[0]))
for ispin in range(nspin):
    for idir in range(4):# 0, Bx, By, Bz
        ndensities[ispin, idir] = numpy.dot(numpy.dot(tmat, densities[ispin, idir]), tmat.transpose())

# ATTENTION values in XDENS for Bx, By, Bz are 2 times bigger for closeshell
# and 4 times biggen for openshell
if nspin == 1:
    ndensities[:, 1:, :, :] *= 2
else:
    ndensities[:, 1:, :, :] *= 4

# open outputfile
outfile=open("XDENS", "w")

# write densities matrices on file
# write in a fortran way
for ispin in range(nspin):
    for idir in range(4):# 0, Bx, By, Bz
        for j in range(ndensities.shape[3]):
            for i in range(ndensities.shape[2]):
                outfile.write("%16.8e\n"%(ndensities[ispin, idir, i, j]))
        outfile.write("\n")
outfile.close()

# write the MOL file with coordinates and basis set
coordinates = None
if hasattr(file, "get_property"):
    coordinates = file.get_property("Current cartesian coordinates").reshape((-1, 3))
nbasisset.write_MOL(filename='MOL', coords=coordinates, turbomole=options.turbomole)


if options.XDENS != 'none':
    if not  os.path.isfile(options.XDENS):
        parser.error("The file ['%s'] does not exist"%(options.XDENS))
    f = open(options.XDENS, 'r')
    lines = f.readlines()
    f.close()
    data = []
    for line in lines:
        if len(line.strip()) != 0:
            data.append(float(line.strip()))
    data = numpy.array(data)
    k = nbasisset.nbasisfunctions
    data = data.reshape((-1, 4, k, k))
    ndensities_o = data.copy()
    # fortran to python arrangement
    for ispin in range(nspin):
        for idir in range(4):
            ndensities_o[ispin, idir] = data[ispin, idir].transpose()
    # check the densities matrices
    for ispin in range(nspin):
        for idir in range(4):
            print("Spin %s, %s"%("alpha + beta" if nspin == 1 else["alpha", "beta"][ispin]\
, "Unperturbed" if idir == 0 else "Magnetic field %i"%(idir)))
            print("        AO orbitals            :   GAUSSIAN     OTHER     DIFF    ")
            for i, ao1 in enumerate(aos.basis):
                for j, ao2 in enumerate(aos.basis):
                    print("%13s -- %13s : %12.8f %12.8f %12.8f"%(ao1.label, ao2.label, ndensities[ispin, idir, j, i], ndensities_o[ispin, idir, j, i], ndensities[ispin, idir, j, i] -ndensities_o[ispin, idir, j, i]))
            print("Maximum Error: %12.8f"%(numpy.amax(numpy.abs(ndensities[ispin, idir] -ndensities_o[ispin, idir]))))
