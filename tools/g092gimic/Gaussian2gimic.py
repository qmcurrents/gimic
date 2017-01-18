#! /usr/bin/env python

import BasisSet
import os.path
from optparse import OptionParser
import numpy
import math
import re

class FCHK (object) :

    def __init__ (self, fchkname='gaussian.fchk'):
        self.filename = fchkname

    def get_property(self,property):
        """
        Read a property name "property" in the fchk file
        return a numpy vector
        """
        f=open(self.filename,'r')
        lines=f.readlines()
        f.close()
        # find the block "property"
        first=0
        last=0
        cols={"R":5,"I":6}
        numpytype={"R":"d","I":"int"}
        stype="R"
        size=0
        for i,l in enumerate(lines):
            if property in l:
                l = l[len(property):]
                if "N=" in l:
                    # read array of data
                    pattern=r"\d+"
                    size=int(re.findall(pattern,l)[0])
                    typepattern=r"   [IR]   "
                    try:
                        stype=re.findall(typepattern,l)[0].strip()
                    except IndexError:
                        raise Exception("The type of data in fchk is not recognized")
                    first=i+1
                    nlines=(size-1)/cols[stype]+1
                    last=first+nlines
                    break
                else:
                    # single data
                    typepattern=r"   [IR]   "
                    try:
                        stype=re.findall(typepattern,l)[0].strip()
                    except IndexError:
                        raise Exception("The type of data in fchk is not recognized")
                    if stype=='I':
                        pattern=r"\d+"
                        return int(re.findall(pattern,l)[0])
                    elif stype=='R':
                        pattern=r"-?\d+\.\d*[ED]?[+\-]?\d*"
                        return float(re.findall(pattern,l)[0])
        lines=lines[first:last]
        if len(lines) ==0:
            return None
        # read the data 
        data=[]
        for line in lines:
            data.extend([float (fl) for fl in line.split()])
        data=numpy.array(data,numpytype[stype])
        if data.size != size:
            raise Exception("The number of data recovered [%i] is not the same as the size written in the file [%i]"%(data.size,size))
        return data

    def get_density(self):
        """
        Get the non-perturbed density
        """
        data = self.get_property("Total SCF Density")
        if type(data) == numpy.ndarray:
            nbasis = self.get_property("Number of basis functions")
            index=0
            density = numpy.zeros((nbasis, nbasis))
            for i in range(nbasis):
                for j in range(i+1):
                    density[i,j]=data[index]
                    index+=1
            
            # fill the other half of the density matrix
            density=density+density.transpose()-numpy.diag(density.diagonal())
        return density


    def get_density_magnetic(self):
        """
        Get the perturbed density by magnetic field
        """
        data = self.get_property("Magnetic Field P1 (GIAO)")
        if type(data) == numpy.ndarray:
            nbasis =  self.get_property("Number of basis functions") 
            index=0
            density = numpy.zeros((3,nbasis, nbasis))
            for idir in range(3):
                for i in range(nbasis):
                    for j in range(i+1):
                        density[idir,i,j]=data[index]
                        index+=1
            
                # fill the other half of the density matrix
                density[idir]=density[idir]-density[idir].transpose()
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
        basisset=[]
        start = 0
        for ishell in range(nshell):
            nprim = primitive_pershell[ishell]
            # get the exponents and coeffs for all the primitive of the shell
            exponents = primitive_exp[start:start+nprim]
            coefficients = primitive_coeff[start:start+nprim]
            coefficients_p =None
            if primitive_coeff_p != None:
                if numpy.nonzero(primitive_coeff_p[start:start+nprim])[0].size !=0 :
                    coefficients_p = primitive_coeff_p[start:start+nprim]
            iatom = shell2atom[ishell]
            an = int(atomicnumbers[iatom-1])
            type = int(shelltype[ishell])
            shell = BasisSet.SHELL(iatom,an,type,exponents=exponents,coefficients=coefficients,coefficients_p=coefficients_p)
            basisset.append(shell)
            start = start + nprim
        if len(basisset) == 0:
            basisset =  None
        else:
            basisset =  BasisSet.BasisSet(basisset)
#            basisset.check_contractedcoeff()
        return basisset



usage = "usage: %prog [options] "
parser = OptionParser(usage=usage)

parser.add_option("-i", "--inputfile", dest="inputfile",
                  help="Input fchk file",
                  action="store", type="string",default="none")

parser.add_option("-t", "--turbomole", dest="turbomole",
                  help="Turbomole XDENS file to compare with the generated one. Default: none",
                  action="store", type="string",default="none")

(options, args) = parser.parse_args()

#check that the user introduce  one input filename
if options.inputfile=="none" :
    parser.error("No inputfile given")

if not os.path.isfile(options.inputfile):
    parser.error("The file ['%s'] does not exist" % options.inputfile)

# open the fchk file 
file=FCHK(options.inputfile)

# read the density matrices
density = file.get_density()
density_m = file.get_density_magnetic()

#read the basis set
basisset = file.get_basisset()

# prepare a transformation matrix from Spherical-harmonic to Cartesian orbitals
# create the block diagonal matrix of transformation
tmat = basisset.Trans2Cart(square=False,order1="gaussian",order2="turbomole",turbomole1=False,turbomole2=True)
tbasisset = basisset.to_Cartesian().split_SP()

# arrange the order of the basis function like turbomole
# S for all atoms, P for all atoms, D for all atoms, ...
aos = tbasisset.to_AOs("turbomole")
index = aos.index_sort_by_type()
if density.shape[0] != tmat.shape[0]:
    raise Exception("The shape of the density matrix [%i] is not the same as the shape of the transformation matrix [%i]"%(density.shape[0],tmat.shape[0]))
if len(index) != tmat.shape[1]:
    raise Exception("The length of the list of index [%i] is not the same as the shape of the transformation matrix [%i]"%(len(index),tmat.shape[1]))
tmat = tmat[:,index]
tbasisset = tbasisset.basis_from_index(tbasisset.index_sort_by_type())
aos = aos.basis_from_index(index)

# apply the transformation to normalized Cartesian 
# with the shells organized by types
# with the angular momenta organized like turbomole
ndensity = numpy.dot(numpy.dot(tmat.transpose(),density),tmat)
ndensity_m = numpy.zeros((3,ndensity.shape[0],ndensity.shape[1]))
for i in range(3):
    ndensity_m[i]  =  numpy.dot(numpy.dot(tmat.transpose(),density_m[i]),tmat)

# ATTENTION values in turbomole are 2 times bigger
ndensity_m *=2

# open outputfile
outfile=open("XDENS","w")

# write density matrix on file
for line in ndensity:
    for val in line:
        outfile.write("%16.8e\n"%(val))
outfile.write("\n")

# write the perturbed density matrices
# write in a fortran way
for idir in range(ndensity_m.shape[0]):
    for j in range(ndensity_m.shape[2]):
        for i in range(ndensity_m.shape[1]):
            outfile.write("%16.8e\n"%(ndensity_m[idir,i,j]))
    outfile.write("\n")

outfile.close()


# write the MOL file with coordinates and basis set
coordinates = file.get_property("Current cartesian coordinates").reshape((-1,3))
tbasisset.write_MOL(filename='MOL',coords=coordinates)

if options.turbomole != 'none':
    if not  os.path.isfile(options.turbomole):
        parser.error("The file ['%s'] does not exist"%(options.turbomole))
    f = open(options.turbomole,'r')
    lines=f.readlines()
    f.close()
    data=[]
    for line in lines:
        if len(line.strip()) !=0:
            data.append(float(line.strip()))
    data = numpy.array(data)
    k = tbasisset.nbasisfunctions
    data = data.reshape((4,k,k))
    ndensity_t = data[0]
    ndensity_m_t = numpy.zeros((3,k,k))
    # fortran to python arrangement
    for i in range(3):
        ndensity_m_t[i] = data[i+1].transpose()
    # check the density matrix
    print "        AO orbitals            :   GAUSSIAN     TURBOMOLE     DIFF    "
    for i,ao1 in enumerate(aos.basis):
        for j,ao2 in enumerate(aos.basis):
            print "%13s -- %13s : %12.8f %12.8f %12.8f"%(ao1.label, ao2.label ,ndensity[j,i],ndensity_t[j,i],ndensity[j,i]-ndensity_t[j,i]) 
    print "Maximum Error: %12.8f"%(numpy.amax(numpy.abs(ndensity-ndensity_t)))

    # check the perturbed density matrix
    for idir in range(3):
        print "Magnetic field %i"%(idir+1)
        print "        AO orbitals            :   GAUSSIAN     TURBOMOLE     DIFF    "
        for i,ao1 in enumerate(aos.basis):
            for j,ao2 in enumerate(aos.basis):
                print "%13s -- %13s : %12.8f %12.8f %12.8f"%(ao1.label,ao2.label,ndensity_m[idir,j,i],ndensity_m_t[idir,j,i],ndensity_m[idir,j,i]-ndensity_m_t[idir,j,i]) 
        print "Maximum Error: %12.8f"%(numpy.amax(numpy.abs(ndensity_m[idir]-ndensity_m_t[idir])))
    
