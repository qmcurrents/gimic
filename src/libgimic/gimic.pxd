cdef extern from "GimicInterface.h":
    cdef cppclass GimicInterface:
        GimicInterface(char *, char *)
        void set_uhf(int)
        void set_magnet(double *b)
        void set_spin(char *s)
        void set_screening(double thrs)
        void calc_jtensor(double *r, double *jt)
        void calc_jvector(double *r, double *jv)

