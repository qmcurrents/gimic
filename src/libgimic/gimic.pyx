from cython.operator cimport dereference as deref
from connector cimport GimicConnector
import numpy as np
cimport gimic

cdef class Gimic(GimicConnector):
    cdef gimic.GimicInterface *thisptr

    def __cinit__(self, mol, xdens):
        if not isinstance(mol, str):
            raise TypeError
        if not isinstance(xdens, str):
            raise TypeError
        self.thisptr = new gimic.GimicInterface(mol, xdens)

    cpdef jtensor(self, r):
        cdef double cr[3]
        cdef double ct[9]
        for i in range(3):
            cr[i] = r[i]
        self.thisptr.calc_jtensor(cr, ct)
        a = np.array(9)
        for i in range(9):
            a[i] = ct[i]
        return a

    cpdef jvector(self, r):
        cdef double cr[3]
        cdef double cv[3]
        for i in range(3):
            cr[i] = r[i]
        self.thisptr.calc_jvector(cr, cv)
        vec=[]
        for i in range(3):
            vec.append(cv[i])
        return vec

    cpdef divj(self, r):
        cdef double cr[3]
        cdef double cd
        for i in range(3):
            cr[i] = r[i]
        self.thisptr.calc_divj(cr, &cd)
        return cd

    cpdef edens(self, r):
        cdef double cr[3]
        cdef double cd
        for i in range(3):
            cr[i] = r[i]
        self.thisptr.calc_edens(cr, &cd)
        return cd

    cpdef set_property(self, prop, val):
        eval('self.set_{0}({1})'.format(prop, repr(val)))

    def set_uhf(self, onoff):
        if not isinstance(onoff, int):
            raise TypeError
        self.thisptr.set_uhf(onoff)

    def set_magnet(self, b):
        cdef double mag[3]
        for i in range(3):
            mag[i] = b[i]
        self.thisptr.set_magnet(mag)

    def set_spin(self, spin):
        if not isinstance(spin, str):
            raise TypeError
        self.thisptr.set_spin(spin)

    def set_screening(self, thrs):
        if not isinstance(thrs, float):
            raise TypeError
        self.thisptr.set_screening(thrs)


