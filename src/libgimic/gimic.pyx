from cython.operator cimport dereference as deref
cimport gimic

cdef class Gimic:
    cdef gimic.GimicInterface *thisptr

    def __cinit__(self, mol, xdens):
        if not isinstance(str, mol):
            raise TypeError
        if not isinstance(str, xdens):
            raise TypeError
        self.thisptr = new gimic.GimicInterface(mol, xdens)

    def set_uhf(self, onoff):
        if not isinstance(int, onoff):
            raise TypeError
        self.thisptr.setUhf(onoff)

    def set_magnet(self, b):
        if not isinstance(list, b):
            raise TypeError
        if not isinstance(float, b[0]):
            raise TypeError
        cdef double mag[3]
        for i in range(3):
            mag[i] = b[i]
        self.thisptr.setMagnet(mag)

    def set_spin(self, spin):
        if not isinstance(str, spin):
            raise TypeError
        self.thisptr.setSpin(spin)

    def set_screening(self, thrs):
        if not isinstance(float, thrs):
            raise TypeError
        self.thisptr.setScreening(thrs)

    def calc_jtensor(self, r):
        if not isinstance(list, r):
            raise TypeError
        if not isinstance(float, r[0]):
            raise TypeError
        cdef double cr[3]
        cdef double ct[9]
        for i in range(3):
            cr[i] = r[i]
        self.thisptr.calcJTensor(cr, ct)
        tens=[]
        for i in range(9):
            tens[i] = ct[i]
        return tens

    def calc_jvector(self, r):
        if not isinstance(list, r):
            raise TypeError
        if not isinstance(float, r[0]):
            raise TypeError
        cdef double cr[3]
        cdef double cv[3]
        for i in range(3):
            cr[i] = r[i]
        self.thisptr.calcJVector(cr, cv)
        vec=[]
        for i in range(3):
            vec[i] = cv[i]
        return vec

    def calc_divj(self, r):
        if not isinstance(list, r):
            raise TypeError
        if not isinstance(float, r[0]):
            raise TypeError
        cdef double cr[3]
        cdef double cd
        for i in range(3):
            cr[i] = r[i]
        self.thisptr.calcDivJ(cr, &cd)
        return cd

    def calc_modj(self, r):
        if not isinstance(list, r):
            raise TypeError
        if not isinstance(float, r[0]):
            raise TypeError
        cdef double cr[3]
        cdef double cd
        for i in range(3):
            cr[i] = r[i]
        self.thisptr.calcModJ(cr, &cd)
        return cd

    def calc_anapole(self, r):
        if not isinstance(list, r):
            raise TypeError
        if not isinstance(float, r[0]):
            raise TypeError
        cdef double cr[3]
        cdef double cd
        for i in range(3):
            cr[i] = r[i]
        self.thisptr.calcAnapole(cr, &cd)
        return cd
