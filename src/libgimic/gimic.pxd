cdef extern from "GimicInterface.h":
    cdef cppclass GimicInterface:
        GimicInterface()
        set_uhf(int)
        set_magnet(double *b)
        set_spin(char *s)
        set_screening(double thrs)
        calc_jtensor(double *r, double *jt)
        calc_jvector(double *r, double *jt)
        calc_divj(double *r, double *dj)
        calc_edens(double *r, double *ed)

