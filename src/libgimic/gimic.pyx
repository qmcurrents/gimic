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
        self.thisptr.set_uhf(onoff)

    def set_magnet(self, b):
        if not isinstance(list, b):
            raise TypeError
        if not isinstance(float, b[0]):
            raise TypeError
        cdef double mag[3]
        for i in range(3):
            mag[i] = b[i]
        self.thisptr.set_magnet(mag)

    def set_spin(self, spin):
        if not isinstance(str, spin):
            raise TypeError
        self.thisptr.set_spin(spin)

    def set_screening(self, thrs):
        if not isinstance(float, thrs):
            raise TypeError
        self.thisptr.set_screening(thrs)

    def calc_jtensor(self, r):
        if not isinstance(list, r):
            raise TypeError
        if not isinstance(float, r[0]):
            raise TypeError
        cdef double cr[3]
        cdef double ct[9]
        for i in range(3):
            cr[i] = r[i]
        self.thisptr.calc_jtensor(cr, ct)
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
        self.thisptr.calc_jvector(cr, cv)
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
        self.thisptr.calc_divj(cr, &cd)
        return cd

    def calc_edens(self, r):
        if not isinstance(list, r):
            raise TypeError
        if not isinstance(float, r[0]):
            raise TypeError
        cdef double cr[3]
        cdef double cd
        for i in range(3):
            cr[i] = r[i]
        self.thisptr.calc_edens(cr, &cd)
        return cd

