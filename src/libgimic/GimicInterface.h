#ifndef GIMIC_INTERFACE_H
#define GIMIC_INTERFACE_H

class GimicInterface {
    public:
        GimicInterface(const char *mol, const char *xdens);
        virtual ~GimicInterface();
        void setUhf(int uhf);
        void setMagnet(const double b[3]);
        void setSpin(char *s);
        void setScreening(double thrs);
        void calcJTensor(const double r[3], double jt[9]);
        void calcJVector(const double r[3], double jv[3]);
        void calcDivJ(const double r[3], double *dj);
        void calcModJ(const double r[3], double *mj);
        void calcAnapole(const double r[3], double *aj);
};
#endif
