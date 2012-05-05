cdef extern from "gausspoints.h": 
    void mkgausspoints(double *a, double *b, int *n, int *o, 
            double *p, double *w)
