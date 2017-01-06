

import numpy
import re
import math

Table=['Xx', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn']
multiplicity={0:1,1:3,-1:4,-2:5,2:6,-3:7,3:10,-4:9,4:15,-5:11,5:21,-6:13,6:28}
mult2type={1:0,3:1,4:-1,5:-2,6:2,7:-3,10:3,9:-4,15:4,11:-5,21:5,13:-6,28:6}
shelltypes = {0:"S",1:"P",2:"D",3:"F",4:"G",5:"H",6:"I",-2:"D",-3:"F",-4:"G",-5:"H",-6:"I",-1:"SP"}
orbitaltypes = ["S","SP","P","D","F","G","H","I"]

angularmomenta_turbomole={}
angularmomenta_turbomole["S"]=[(0,0,0)]
angularmomenta_turbomole["P"]=[(1,0,0),(0,1,0),(0,0,1)]
angularmomenta_turbomole["SP"]=[(0,0,0),(1,0,0),(0,1,0),(0,0,1)]
angularmomenta_turbomole["D"]=[(2,0,0),(0,2,0),(0,0,2),(1,1,0),(1,0,1),(0,1,1)]
angularmomenta_turbomole["F"]=[(3,0,0),(0,3,0),(0,0,3),(2,1,0),(2,0,1),(1,2,0),(0,2,1),(1,0,2),(0,1,2),(1,1,1)]
angularmomenta_turbomole["G"]=[(4,0,0),(0,4,0),(0,0,4),(3,1,0),(3,0,1),(1,3,0),(0,3,1),(1,0,3),(0,1,3),(2,2,0),(2,0,2),(0,2,2),(2,1,1),(1,2,1),(1,1,2)]
angularmomenta_turbomole["H"]=[(5,0,0),(0,5,0),(0,0,5),(4,1,0),(4,0,1),(1,4,0),(0,4,1),(1,0,4),(0,1,4),(3,2,0),(3,0,2),(2,3,0),(0,3,2),(2,0,3),(0,2,3),(3,1,1),(1,3,1),(1,1,3),(2,2,1),(2,1,2),(1,2,2)]
angularmomenta_turbomole["I"]=[(6,0,0),(0,6,0),(0,0,6),(5,1,0),(5,0,1),(1,5,0),(0,5,1),(1,0,5),(0,1,5),(4,2,0),(4,0,2),(2,4,0),(0,4,2),(2,0,4),(0,2,4),(4,1,1),(1,4,1),(1,1,4),(3,3,0),(3,0,3),(0,3,3),(3,2,1),(3,1,2),(2,3,1),(1,3,2),(2,1,3),(1,2,3),(2,2,2)]

angularmomenta_gaussian={}
angularmomenta_gaussian["S"]=[(0,0,0)]
angularmomenta_gaussian["P"]=[(1,0,0),(0,1,0),(0,0,1)]
angularmomenta_gaussian["SP"]=[(0,0,0),(1,0,0),(0,1,0),(0,0,1)]
angularmomenta_gaussian["D"]=[(2,0,0),(0,2,0),(0,0,2),(1,1,0),(1,0,1),(0,1,1)]
angularmomenta_gaussian["F"]=[(3,0,0),(0,3,0),(0,0,3),(1,2,0),(2,1,0),(2,0,1),(1,0,2),(0,1,2),(0,2,1),(1,1,1)]
angularmomenta_gaussian["G"]=[(0,0,4),(0,1,3),(0,2,2),(0,3,1),(0,4,0),(1,0,3),(1,1,2),(1,2,1),(1,3,0),(2,0,2),(2,1,1),(2,2,0),(3,0,1),(3,1,0),(4,0,0)]
angularmomenta_gaussian["H"]=[(0,0,5),(0,1,4),(0,2,3),(0,3,2),(0,4,1),(0,5,0),(1,0,4),(1,1,3),(1,2,2),(1,3,1),(1,4,0),(2,0,3),(2,1,2),(2,2,1),(2,3,0),(3,0,2),(3,1,1),(3,2,0),(4,0,1),(4,1,0),(5,0,0)]
angularmomenta_gaussian["I"]=[(0,0,6),(0,1,5),(0,2,4),(0,3,3),(0,4,2),(0,5,1),(0,6,0),(1,0,5),(1,1,4),(1,2,3),(1,3,2),(1,4,1),(1,5,0),(2,0,4),(2,1,3),(2,2,2),(2,3,1),(2,4,0),(3,0,3),(3,1,2),(3,2,1),(3,3,0),(4,0,2),(4,1,1),(4,2,0),(5,0,1),(5,1,0),(6,0,0)]


Fields = ["x","y","z"]

mqn={}
mqn["S"]=[0]
mqn["SP"]=[0,1,-1,0]
mqn["P"]=[1,-1,0]
mqn["D"]=[0,1,-1,2,-2]
mqn["F"]=[0,1,-1,2,-2,3,-3]
mqn["G"]=[0,1,-1,2,-2,3,-3,4,-4]
mqn["H"]=[0,1,-1,2,-2,3,-3,4,-4,5,-5]
mqn["I"]=[0,1,-1,2,-2,3,-3,4,-4,5,-5,6,-6]


def coeff_Cart2Spher(m,l,real=False,turbomole=False):
    """
    Coefficient of transformation between a normalized Cartesian orbitals and a normalized spherical-harmonic ones
    \chi^sphe_{l,m} = \sum_{lx,ly,lz} c_{l,m,lx,ly,lz) \chi^cart_{lx,ly,lz}
    Param m: m value for normalized spherical-harmonic orbital
    Param l: l vector with lx,ly,lz values for normalized Cartesian orbital
    Param real: transform into real spherical-harmonic instead of complex (default)
    Param turbomole: use turbomole normalization for Cartesian orbitals
    """
    lx,ly,lz = l
    l = sum(l)
    j = (lx + ly - abs(m)) 
    if j%2 == 1:
        return 0.0
    j = j // 2
    N = math.sqrt( float(math.factorial(2*lx) * math.factorial(2*ly) * math.factorial(2*lz) * math.factorial(l) * math.factorial(l-abs(m)) ) /\
                 float( math.factorial(2*l) * math.factorial(lx) * math.factorial(ly) * math.factorial(lz) * math.factorial(l+abs(m)) ) )
    N = N / (2**l * math.factorial(l))
    fact1 = 0.0
    for i in range((l-abs(m))//2+1):
        bin1 = math.factorial(l) / ( math.factorial(i) * math.factorial(l - i))
        try:
            bin2 = math.factorial(i) / ( math.factorial(j) * math.factorial(i - j))
        except ValueError:
            bin2 = 0.0
        N2 = (-1)**i * math.factorial(2*l-2*i) / math.factorial(l - abs(m) - 2*i)
        fact1 = fact1 + bin1*bin2*N2
    fact2 = 0.0
    for k in range(j+1):
        bin3 = math.factorial(j) / ( math.factorial(k) * math.factorial(j - k))
        try:
            bin4 = math.factorial(abs(m)) / ( math.factorial(lx-2*k) * math.factorial(abs(m) - lx + 2*k)) 
        except ValueError:
            bin4 = 0.0
        N3 = (abs(m)-lx+2*k)
#        if N3 % 2 == 0:
#            N3 = (-1)**(N3//2)
#        else:
#            N3 = numpy.sign(m)*1j *( (-1)**((N3-1)//2))
        if not real:
            # N3 = (-1)**(\pm N3/2)
            fp = [1,1j,-1,-1j]
            fm = [1,-1j,-1,1j]
        else:
            # N3 = [ (-1)**(N3/2) + (-1)**(-N3/2)] / \sqrt{2}  for m>0
            # N3 = [ (-1)**(N3/2) - (-1)**(-N3/2)] / \sqrt{2}i for m<0
            fp = [math.sqrt(2.0),0,-math.sqrt(2.0),0]
            fm = [0,math.sqrt(2.0),0,-math.sqrt(2.0)]
        if m > 0:
            N3 = fp[N3 % 4]
        elif m < 0:
            N3 = fm[N3 % 4]
        else:
            N3 = 1.0
        fact2 = fact2 + bin3 * bin4 * N3
    coeff =  N * fact1 * fact2
    if turbomole:
        coeff = coeff * math.sqrt(float( (math.factorial(lx)*math.factorial(ly)*math.factorial(lz))**3 * 2**l ) / float( math.factorial(2*lx) * math.factorial(2*ly) * math.factorial(2*lz) ))
    return coeff


class SHELL(object):
    """
    Class to describe one shell 
    A shell is described by a contraction of GTO
    A shell is quantified by l quantum number
    """

    def __init__(self,iatom,atnum,l,exponents,coefficients,coefficients_p=None,coord=None):
        """
        Param iatom: atom number (one-based)
        Param atnum: atomic number 
        Param l: int to represent the type of atomic orbitals
                    the integer correspond to the 'l' quantum number with a sign
                    + sign means Cartesian function
                    - sign means Spherical-harmonic function
                    -1 means SP orbital
        Param exponents: numpy array with exponents of the contraction of GTO
        Param coefficients: numpy array with coefficients of the contraction of GTO
        Param coefficients_p: numpy array with coefficients of p orbitals for SP type 
        Param coord: numpy array with the coordinates of the atoms in Angstrom
        """
        self.iatom = iatom
        self.atnum = atnum
        self.type = l
        self.exponents = exponents
        self.coefficients = coefficients
        self.coefficients_p = coefficients_p
        self.coord = coord
        self.nprimitive = self.exponents.size
        self.multiplicity = multiplicity[self.type]
        self.shelltype = shelltypes[self.type]
        self.label = Table[self.atnum]+"%i@"%(self.iatom)+self.shelltype

    def copy(self):
        newSHELL=SHELL(self.iatom,self.atnum,self.type,self.exponents,self.coefficients,self.coefficients_p,self.coord)
        return newSHELL

    def to_Cartesian(self):
        """
        Return a new SHELL where the orbitals are Cartesian
        """
        type = self.type
        if type < -1:
            type = -type
        newSHELL=SHELL(self.iatom,self.atnum,type,self.exponents,self.coefficients,self.coefficients_p,self.coord)
        return newSHELL

    def to_Spherical(self):
        """
        Return a new SHELL where the orbitals are Spherical-harmonic
        """
        type = self.type
        if type > 1:
            type = -type
        newSHELL=SHELL(self.iatom,self.atnum,type,self.exponents,self.coefficients,self.coefficients_p,self.coord)
        return newSHELL


    @staticmethod
    def mat_Cart2Spher(shelltype,square=False,real=False,order="gaussian",turbomole=False):
        """
        Transformation matrices from normalized Cartesian orbitals to normalized spherical-harmonic ones
        for one type of orbital, one shell
        \chi^sphe_{l,m} = \sum_{lx,ly,lz} c_{l,m,lx,ly,lz) \chi^cart_{lx,ly,lz}
        Lines : normalized cartesian orbitals. Columns: normalized spherical-harmonic orbitals
        Param shelltype: one of orbitaltypes
        Param square: return square matrices with blank lines instead of rectangular matrices (default)
        Param real: transform into real spherical-harmonic instead of complex (default)
        Param order: order of the angularmomenta. Either "gaussian" or "turbomole"
        Param turbomole: use turbomole normalization for Cartesian orbitals
        """
        mqns = mqn[shelltype]
        angmoms = eval("angularmomenta_"+order)[shelltype]
        if real:
            dtype = "float"
        else:
            dtype = "complex"
        if square:
            c = numpy.zeros((len(angmoms),len(angmoms)),dtype)
        else:
            c = numpy.zeros((len(angmoms),len(mqns)),dtype)
        for il,l in enumerate(angmoms):
            for im,m in enumerate(mqns):
                c[il,im]= coeff_Cart2Spher(m,l,real=real,turbomole=turbomole)
        return c

    @staticmethod
    def Overlap_Cart(shelltype,order="gaussian",turbomole=False):
        """
        Overlap matrix between Cartesian atomic orbitals of same l
        Param shelltype: one of orbitaltypes
        Param order: order of the angularmomenta. Either "gaussian" or "turbomole"
        Param turbomole: use turbomole normalization for Cartesian orbitals
        """
        angmoms = eval("angularmomenta_"+order)[shelltype]
        ov = numpy.zeros((len(angmoms),len(angmoms)))
        for il,l1 in enumerate(angmoms):
            for jl,l2 in enumerate(angmoms):
                if len(numpy.flatnonzero( (numpy.array(l1)+numpy.array(l2)) % 2 )) >0:
                    ov[il,jl] = 0.0
                else:
                    num1 = math.factorial(l1[0]+l2[0]) * math.factorial(l1[1]+l2[1]) * math.factorial(l1[2]+l2[2])
                    den1 = math.factorial((l1[0]+l2[0])//2) * math.factorial((l1[1]+l2[1])//2) * math.factorial((l1[2]+l2[2])//2)
                    num2 = math.factorial(l1[0]) * math.factorial(l1[1]) * math.factorial(l1[2]) *\
                        math.factorial(l2[0]) * math.factorial(l2[1]) * math.factorial(l2[2]) 
                    den2 = math.factorial(2*l1[0]) * math.factorial(2*l1[1]) * math.factorial(2*l1[2]) *\
                        math.factorial(2*l2[0]) * math.factorial(2*l2[1]) * math.factorial(2*l2[2]) 
                    ov[il,jl] = float(num1) / float(den1) * math.sqrt(float(num2)/float(den2))
                    if turbomole:
                        f1 = math.sqrt( float(math.factorial(2*l1[0]) * math.factorial(2*l1[1]) * math.factorial(2*l1[2]) ) / float( (math.factorial(l1[0])*math.factorial(l1[1])*math.factorial(l1[2]))**3 * 2**(sum(l1)) ))
                        f2 = math.sqrt( float(math.factorial(2*l2[0]) * math.factorial(2*l2[1]) * math.factorial(2*l2[2]) ) / float( (math.factorial(l2[0])*math.factorial(l2[1])*math.factorial(l2[2]))**3 * 2**(sum(l2)) ))
                        ov[il,jl] = ov[il,jl] * f1 * f2
        return ov

    @staticmethod
    def mat_Spher2Cart(shelltype,square=False,real=False,order="gaussian",turbomole=False):
        """
        Transformation matrices from normalized spherical-harmonic orbitals to normalized Cartesian ones
        for one type of orbital
        \chi^cart_{lx,ly,lz} = \sum_{l<=lx+ly+lz} c-1_{l,m,ly,lz,lz) \chi^sphe_{l,m}
        c-1 = c^dagger S since c^dagger S c = 1 = c-1 c
        Lines: normalized spherical-harmonic orbitals. Columns : normalized cartesian orbitals. 
        Param shelltype: one of orbitaltypes
        Param square: return square matrices with blank lines instead of rectangular matrices (default)
        Param real: transform into real spherical-harmonic instead of complex (default)
        Param order: order of the angularmomenta. Either "gaussian" or "turbomole"
        Param turbomole: use turbomole normalization for Cartesian orbitals
        """
        C = SHELL.mat_Cart2Spher(shelltype,square=square,real=real,order=order,turbomole=turbomole)
        if numpy.dtype("complex") == C.dtype:
            C = C.conjugate()
        S = SHELL.Overlap_Cart(shelltype,order=order,turbomole=turbomole)
        Cm1 = numpy.dot(C.transpose(),S)
        return Cm1

    @staticmethod
    def mat_Cart2Cart(shelltype,order1="gaussian",order2="gaussian",turbomole1=False,turbomole2=False):
        """
        Transformation matrices from normalized Cartesian orbitals to normalized Cartesian ones
        for one type of orbital
        The arrangement of the angular momenta of the orbitals can be different between the two sets
        Param shelltype: one of orbitaltypes
        Param order1, order1: order of the angularmomenta. Either "gaussian" or "turbomole"
        Param turbomole1, turbomole2: use turbomole normalization for Cartesian orbitals
        """
        angmoms1 = eval("angularmomenta_"+order1)[shelltype]
        angmoms2 = eval("angularmomenta_"+order2)[shelltype]
        c = numpy.zeros((len(angmoms1),len(angmoms2)))
        for il,l1 in enumerate(angmoms1):
            fct = 1.0
            if turbomole1 and (not turbomole2):
                fct = math.sqrt( float( (math.factorial(l1[0])*math.factorial(l1[1])*math.factorial(l1[2]))**3 * 2**(sum(l1)) ) / float(math.factorial(2*l1[0]) * math.factorial(2*l1[1]) * math.factorial(2*l1[2]) ) )
            elif (not turbomole1) and turbomole2:
                fct = math.sqrt( float(math.factorial(2*l1[0]) * math.factorial(2*l1[1]) * math.factorial(2*l1[2]) ) / float( (math.factorial(l1[0])*math.factorial(l1[1])*math.factorial(l1[2]))**3 * 2**(sum(l1)) ))
            jl = angmoms2.index(l1)
            c[il,jl]= fct 
        return c

    @staticmethod
    def mat_Spher2Spher(shelltype,real1=False,real2=False):
        """
        Transformation matrices from normalized spherical-harmonic orbitals to normalized spherical-harmonic ones
        for one type of orbital
        Both sets can be real or complex
        Param shelltype: one of orbitaltypes
        Param real1, real1: transform into real spherical-harmonic instead of complex (default)
        """
        mqns = mqn[shelltype]
        if (real1 == real2):
            c = numpy.eye(len(mqns),len(mqns))
        else:
            # transformation matrix from imaginary to real
            ci2r = numpy.zeros((len(mqns),len(mqns)),"complex")
            for im,m1 in enumerate(mqns):
                if m1 == 0:
                    ci2r[im,im]= 1.0
                elif m1 > 0:
                    ci2r[im,im]= 1.0/math.sqrt(2.0)
                    jm = mqns.index(-m1)
                    ci2r[jm,im]= 1.0/math.sqrt(2.0)
                else:
                    ci2r[im,im]= 1.0j/math.sqrt(2.0)
                    jm = mqns.index(-m1)
                    ci2r[jm,im]= -1.0j/math.sqrt(2.0)
            if real2:
                c = ci2r
            else:
                c = ci2r.conjugate().transpose()
        return c

    def split_SP(self):
        """
        Split a SP shell in a S and a P shells
        return a list of two shells out of the SP shell 
               or a list with one shell as a copy
        """
        shells = []
        type=self.type
        if type == -1:
            shells.append(SHELL(self.iatom,self.atnum,0,exponents=self.exponents,coefficients=self.coefficients,coefficients_p=None,coord=self.coord))
            shells.append(SHELL(self.iatom,self.atnum,1,exponents=self.exponents,coefficients=self.coefficients_p,coefficients_p=None,coord=self.coord))
        else:
            shells.append(self.copy())
        return shells

    def to_AO(self,order):
        """
        Convert a shell to a list of AOs
        Param order: order of the angularmomenta. Either "gaussian" or "turbomole"
        """
        AOs = []
        type = self.type
        if type == -1: #special for SP
            AOs.append(AO(self.iatom,self.atnum,-1,(0,0,0),exponents=self.exponents,coefficients=self.coefficients,coord=self.coord))
            AOs.append(AO(self.iatom,self.atnum,-1,(1,0,0),exponents=self.exponents,coefficients=self.coefficients_p,coord=self.coord))
            AOs.append(AO(self.iatom,self.atnum,-1,(0,1,0),exponents=self.exponents,coefficients=self.coefficients_p,coord=self.coord))
            AOs.append(AO(self.iatom,self.atnum,-1,(0,0,1),exponents=self.exponents,coefficients=self.coefficients_p,coord=self.coord))
        else:
            for im in range(self.multiplicity):
                if type < -1:
                    m= mqn[self.shelltype][im]
                else:
                    m = eval("angularmomenta_"+order)[self.shelltype][im]
                AOs.append(AO(self.iatom,self.atnum,type,m,exponents=self.exponents,coefficients=self.coefficients,coord=self.coord))
        return AOs


class AO(object):
    """
    Class to describe a AO function
    An AO function is described by a contraction of GTO
    An AO is quantified by l and m quantum numbers
    """

    def __init__(self,iatom,atnum,l,m,exponents,coefficients,coord=None):
        """
        Param iatom: atom number (one-based)
        Param atnum: atomic number 
        Param l: int to represent the type of atomic orbitals
                    the integer correspond to the 'l' quantum number with a sign
                    + sign means Cartesian function
                    - sign means Spherical-harmonic function
                    -1 means SP orbital
        Param m: m quantum number either a int (spherical-harmonic) or tuple (cartesian)
        Param exponents: numpy array with exponents of the contraction of GTO
        Param coefficients: numpy array with coefficients of the contraction of GTO
        Param coord: numpy array with the coordinates of the atoms in Angstrom
        """
        self.iatom = iatom
        self.atnum = atnum
        self.type = l
        self.m = m
        self.exponents = exponents
        self.coefficients = coefficients
        self.coefficients_p = None
        self.coord = coord
        self.nprimitive = self.exponents.size
        self.multiplicity = 1
        self.shelltype = shelltypes[self.type]
        self.label = Table[self.atnum]+"%i@"%(self.iatom)+self.shelltype
        if type(self.m) == int:
            self.label = self.label + "%i"%(self.m)
        elif len(self.m)== 3:
            power=""
            for idir in range(3):
                for i in range(self.m[idir]):
                    power += Fields[idir]
            self.label = self.label + power
        else:
            raise Exception("m quantum number for AO is not correct")
            

    def copy(self):
        newAO=AO(self.iatom,self.atnum,self.type,self.m,self.exponents,self.coefficients,self.coord)
        return newAO

class BasisSet(object):

    def __init__(self,SHELLs):
        """
        Param SHELLs: list of SHELL objects
        """
        
        self.basis = SHELLs
        self.nshells=len(self.basis)
        self.maptype=[]
        self.mapatom=[]
        self.mapatnum=[]
        self.nbasisfunctions=0
        for shell in self.basis:
            self.maptype.append(shell.type)
            self.mapatom.append(shell.iatom)
            self.mapatnum.append(shell.atnum)
            self.nbasisfunctions +=shell.multiplicity
        self.maptype = numpy.array(self.maptype)
        self.mapatom = numpy.array(self.mapatom)
        self.mapatnum = numpy.array(self.mapatnum)
        self.atoms   = sorted(set(self.mapatom)) # list of atoms
        self.atnums  = sorted(set(self.mapatnum)) #list of atnums
        self.natoms = len(self.atoms)

    def to_Cartesian(self):
        """
        Return a new basis set where all the orbitals are Cartesian
        """
        newbasis=[shell.to_Cartesian() for shell in self.basis]
        newbasis = BasisSet(newbasis)
        return newbasis

    def to_Spherical(self):
        """
        Return a new basis set where all the orbitals are Spherical-harmonic
        """
        newbasis=[shell.to_Spherical() for shell in self.basis]
        newbasis = BasisSet(newbasis)
        return newbasis

    def split_SP(self):
        """
        Return a new basis set where the SP functions have been converted into S and P shells
        """
        if not isinstance(self.basis[0],SHELL):
            raise Exception("This method is only available for basis set made of shells")
        newbasis=[]
        for shell in self.basis:
            newbasis.extend(shell.split_SP())
        newbasis = BasisSet(newbasis)
        return newbasis

    def to_AOs(self,order):
        """
        Return a new basis set where the list is made of AO functions instead of shells
        Param order: order of the angularmomenta. Either "gaussian" or "turbomole"
        """
        if isinstance(self.basis[0],SHELL):
            newbasis=[]
            for shell in self.basis:
                newbasis.extend(shell.to_AO(order=order))
        else:
            newbasis=[ ao.copy() for ao in self.basis]
        newbasis = BasisSet(newbasis)
        return newbasis

    def Trans2Spher(self,square=False,order="gaussian",turbomole=False):
        """
        Transformation matrices to normalized real spherical-harmonic orbitals for a basis set
        Lines : normalized orbitals. Columns: normalized spherical-harmonic orbitals
        Param square: return square matrices with blank lines instead of rectangular matrices (default)
        Param order: order of the angularmomenta. Either "gaussian" or "turbomole"
        Param turbomole: use turbomole normalization for Cartesian orbitals
        """
        index1 = 0
        index2 = 0
        m = numpy.zeros((10000,10000))
        # Build the block-diagonal transformation matrix m
        for shell in self.basis:
            if shell.type <= 1:
                # already in spherical
                data = SHELL.mat_Spher2Spher(shell.shelltype,real1=True,real2=True)
            else:
                data = SHELL.mat_Cart2Spher(shell.shelltype,square=square,real=True,order=order,turbomole=turbomole)
            end1 = index1 + data.shape[0]
            end2 = index2 + data.shape[1]
            m[index1:end1,index2:end2] = data
            index1 = end1
            index2 = end2
        return m[:index1,:index2]

    def Trans2Cart(self,square=False,order1="gaussian",order2="gaussian",turbomole1=False, turbomole2=False):
        """
        Transformation matrices to normalized Cartesian orbitals for a basis set
        Lines: normalized orbitals. Columns : normalized cartesian orbitals. 
        Param square: return square matrices with blank lines instead of rectangular matrices (default)
        Param order: order of the angularmomenta. Either "gaussian" or "turbomole"
        Param turbomole: use turbomole normalization for Cartesian orbitals
        """
        index1 = 0
        index2 = 0
        m = numpy.zeros((10000,10000))
        # Build the block-diagonal transformation matrix m
        for shell in self.basis:
            if shell.type >= -1:
                # already in Cartesian
                data = SHELL.mat_Cart2Cart(shell.shelltype,order1=order1,order2=order2,turbomole1=turbomole1,turbomole2=turbomole2)
            else:
                data = SHELL.mat_Spher2Cart(shell.shelltype,square=square,real=True,order=order2,turbomole=turbomole2)
            end1 = index1 + data.shape[0]
            end2 = index2 + data.shape[1]
            m[index1:end1,index2:end2] = data
            index1 = end1
            index2 = end2
        return m[:index1,:index2]

    def get_index_for_atom(self,iatom):
        """
        return the index of the SHELLs associated to an atom 'iatom'
        Param iatom: atom number (one-based)
        """
        index = numpy.where(self.mapatom == iatom)[0]
        return index

    def get_index_for_type(self,type):
        """
        return the index of the SHELLs associated to a type 
        Param type: int to represent the type of atomic orbitals
                    the integer correspond to the 'l' quantum number with a sign
                    + sign means Cartesian function
                    - sign means Spherical-harmonic function
                    -1 means SP orbital
        """
        index = numpy.where(self.maptype == type)[0]
        return index

    def index_sort_by_type(self):
        """
        Sort the SHELLs by type (slow axis) and by atom (fast axis)
        """
        mapatom = self.mapatom
        # convert self.maptype into index array using 'orbitaltypes' order
        maptype = numpy.array([orbitaltypes.index(shelltypes[t]) for t in self.maptype])
        # first sort along mapatom
        index1 = numpy.argsort(mapatom,kind="mergesort")
        # sort maptype with these index
        maptype = maptype[index1]
        # then sort along maptype
        index2 = numpy.argsort(maptype,kind="mergesort")
        index= index1[index2]
        return index

    def index_sort_by_atom(self):
        """
        Sort the SHELLs by atom (slow axis) and by type (fast axis)
        """
        mapatom = self.mapatom
        # convert self.maptype into index array using 'orbitaltypes' order
        maptype = numpy.array([orbitaltypes.index(shelltypes[t]) for t in self.maptype])
        # first sort along maptype
        index1 = numpy.argsort(maptype,kind="mergesort")
        # sort mapatom with these index
        mapatom = mapatom[index1]
        # then sort along mapatom
        index2 = numpy.argsort(mapatom,kind="mergesort")
        index= index1[index2]
        return index

    def basis_from_index(self,index):
        """
        Return a new basis set with the SHELLs in index
        Param index: list of index
        """
        basis = BasisSet([ self.basis[i].copy() for i in index])
        return basis


    def get_atombasis(self,iatom):
        """
        Return a new basis set with the basis functions associated to an atom iatom
        """
        return self.basis_from_index(self.get_index_for_atom(iatom))

    def get_nborbitals(self):
        """
        Return a list with the number orbitals for each 'orbitaltypes'
        """
        typelist = [[k for k,v in shelltypes.items() if v==t] for t in orbitaltypes]
        types=[]
        for ts in typelist:
            index=[]
            for t in ts:
                index.extend(self.get_index_for_type(t))
            types.append(len(index))
        return types

    def write(self,file):
        """
        write the basis set in file
        """
        for atnum in self.atnums:
            iatom = self.mapatom[numpy.where(self.mapatnum == atnum)][0]
            atombasis = self.get_atombasis(iatom)
            file.write("****\n")
            file.write("%s\n"%(Table[atnum]))
            file.write("Nb_shells %i\n"%(atombasis.nshells))
            for shell in atombasis.basis:
                file.write("%s  %i \n"%(shell.shelltype,shell.nprimitive))
                coeffs = shell.coefficients
                coeffs_p = shell.coefficients_p
                exp = shell.exponents
                if coeffs_p!= None:
                    for a,b,c in zip(exp,coeffs,coeffs_p):
                        file.write("  %18.10E  %18.10E  %18.10E\n"%(a,b,c))
                else:
                     for a,b in zip(exp,coeffs):
                        file.write("  %18.10E  %18.10E\n"%(a,b))
                    
    def write_MOL(self,filename,coords=None):
        """
        write the basis set in file
        Param filename: filename to write into
        Param coords: coordinates of the atoms
        """
        HEADER="""INTGRL        1    0    1    0    0    0    0    0    0
TURBOMOLE
              Generated by Gaussian2gimic
"""
        file=open(filename,"w")
        file.write(HEADER)
        file.write("%i    0            0.10E-08              0    0\n"%(self.natoms))
        file.write("9999.00      3.00\n")
        for iatom  in self.atoms:
            atombasis=self.get_atombasis(iatom)
            nborbitals = [x for x in atombasis.get_nborbitals() if x != 0] # keep the nb of orbitals that are non-zero
            atnum = atombasis.atnums[0]
            if coords != None:
                coord = coords[iatom-1]
            else:
                coord = atombasis.basis[0].coord
                if coord == None:
                    raise Exception("No atomic coordinate found for atom %i"%(iatom))
            file.write("%2.1f  %3i %3i"%(float(atnum),1,len(nborbitals)))
            for norb in nborbitals:
                file.write(" %3i"%(norb))
            file.write("\n")
            file.write("%s %i %19.12f %19.12f %19.12f\n"%(Table[atnum],1,coord[0],coord[1],coord[2]))
            for shell in atombasis.basis:
                file.write("   %3i %3i\n"%(shell.nprimitive,1))
                coeffs = shell.coefficients
                coeffs_p = shell.coefficients_p
                exp = shell.exponents
                if coeffs_p!= None:
                    for a,b,c in zip(exp,coeffs,coeffs_p):
                        file.write("  %16.10f  %16.10f  %16.10f\n"%(a,b,c))
                else:
                     for a,b in zip(exp,coeffs):
                        file.write("  %16.10f  %16.10f\n"%(a,b))
        file.close()

    @staticmethod
    def read_MOL(filename):
        """
        read the basis set from MOL file
        Param filename: filename of the MOL file
        """
        file=open(filename,"r")
        lines = file.readlines()
        file.close()
        natoms = int(lines[3].split()[0])
        basisset=[]
        lines = lines[5:]
        iatom=1
        # remove empty lines
        lines = [ l for l in lines if len(l.strip())!=0]
        while(len(lines)>0):
            data = lines[0].strip().split()
            atnum = int(float(data[0]))
            nborbitals=[int(x) for x in data[3:]]
            coord= [float(x) for x in lines[1].split()[2:]]
            lines = lines[2:]
            for itype,nborbital in enumerate(nborbitals):
                for iorbital in range(nborbital):
                    nprim = int(lines[0].split()[0])
                    block = lines[1:1+nprim]
                    lines = lines[1+nprim:]
                    data = []
                    for line in block:
                        data.append([float(x) for x in line.split()])
                    data = numpy.array(data)
                    exponents=data[:,0]
                    coefficients=data[:,1]
                    try:
                        coefficients_p=data[:,2]
                    except IndexError:
                        coefficients_p=None
                    # FIXME: itype does not say if Cart or Spherical and no SP
                    shell = SHELL(iatom,atnum,itype,exponents,coefficients,coefficients_p,coord=coord)
                    basisset.append(shell)
            iatom +=1
        if len(basisset) == 0:
            basisset =  None
        else:
            basisset =  BasisSet(basisset)
        return basisset

    def check_contractedcoeff(self):
        """
        Check that=
        \sum_p,q d_p\mu d_q\mu  2^(l+3/2)  \sqrt(a_p^l+3/2 a_q^l+3/2) / (a_p+a_q)^l+3/2
        """
        for shell in self.basis:
            type = shell.type
            if type == -1:
                type=[0,1]
            else:
                type=[abs(type)]
            coeffs = shell.coefficients
            coeffs_p = shell.coefficients_p
            exponents = shell.exponents
            if coeffs_p == None:
                coefficients = coeffs[numpy.newaxis,:]
            else:
                coefficients = numpy.vstack((coeffs,coeffs_p))
            # loop over all set of coefficients for each l 
            # (only important for sp shell)
            value=[]
            for il,l in enumerate(type):
                aa = numpy.sqrt(numpy.outer(exponents,exponents)**(l+1.5))
                dd = numpy.outer(coefficients[il],coefficients[il])
                apa = numpy.add.outer(exponents,exponents)**(l+1.5)
                value.append( (2**(l+1.5) * dd * aa / apa).sum() )
            print shell.label,["%.8f"%(v) for v in value]

