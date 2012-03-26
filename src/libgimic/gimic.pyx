from cython.operator cimport dereference as deref
cimport gimic

cdef class Gimic:
    cdef gimic.GimicInterface *thisptr

    def __cinit__(self, mol, xdens):
        if not isinstance(mol, str):
            raise TypeError
        if not isinstance(xdens, str):
            raise TypeError
        self.thisptr = new gimic.GimicInterface(mol, xdens)

    def set_uhf(self, onoff):
        if not isinstance(onoff, int):
            raise TypeError
        self.thisptr.set_uhf(onoff)

    def set_magnet(self, b):
        if not isinstance(b, list):
            raise TypeError
        if not isinstance(b[0], float):
            raise TypeError
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

    def calc_jtensor(self, r):
        if not isinstance(r, list):
            raise TypeError
        if not isinstance(r[0], float):
            raise TypeError
        cdef double cr[3]
        cdef double ct[9]
        for i in range(3):
            cr[i] = r[i]
        self.thisptr.calc_jtensor(cr, ct)
        tens=[]
        for i in range(9):
            tens.append(ct[i])
        return tens

    def calc_jvector(self, r):
        if not isinstance(r, list):
            raise TypeError
        if not isinstance(r[0], float):
            raise TypeError
        cdef double cr[3]
        cdef double cv[3]
        for i in range(3):
            cr[i] = r[i]
        self.thisptr.calc_jvector(cr, cv)
        vec=[]
        for i in range(3):
            vec.append(cv[i])
        return vec

    def calc_divj(self, r):
        if not isinstance(r, list):
            raise TypeError
        if not isinstance(r[0], float):
            raise TypeError
        cdef double cr[3]
        cdef double cd
        for i in range(3):
            cr[i] = r[i]
        self.thisptr.calc_divj(cr, &cd)
        return cd

    def calc_edens(self, r):
        if not isinstance(r, list):
            raise TypeError
        if not isinstance(r[0], float):
            raise TypeError
        cdef double cr[3]
        cdef double cd
        for i in range(3):
            cr[i] = r[i]
        self.thisptr.calc_edens(cr, &cd)
        return cd

