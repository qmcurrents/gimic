cdef extern from "gausspoints.h": 
    void mkgausspoints(double *a, double *b, int *n, double *p, double *w)
