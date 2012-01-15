cdef extern from "GimicInterface.h":
    cdef cppclass GimicInterface:
        GimicInterface(char *, char *)
        void setUhf(int)
        void setMagnet(double *b)
        void setSpin(char *s)
        void setScreening(double thrs)
        void calcJTensor(double *r, double *jt)
        void calcJVector(double *r, double *jt)
        void calcDivJ(double *r, double *dj)
        void calcModJ(double *r, double *mj)
        void calcAnapole(double *r, double *aj)

