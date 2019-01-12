#ifndef GIMIC_INTERFACE_H
#define GIMIC_INTERFACE_H

class GimicInterface {
    public:
        GimicInterface(const char *mol, const char *xdens);
        virtual ~GimicInterface();
        void set_uhf(int uhf);
        void set_magnet(const double b[3]);
        void set_spin(char *s);
        void set_screening(double thrs);
        void calc_jtensor(const double r[3], double jt[9]);
        void calc_jvector(const double r[3], double jv[3]);
        void calc_modj(const double r[3], double *mj);
};
#endif
