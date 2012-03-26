#include <iostream>
#include "GimicInterface.h"
#include "gimic_interface.h"

using namespace std;

GimicInterface::GimicInterface(const char *mol, const char *xdens) {
    gimic_init(mol, xdens);
}

GimicInterface::~GimicInterface() {
    gimic_finalize();
}

void GimicInterface::setUhf(int uhf) {
    gimic_set_uhf(&uhf);
}

void GimicInterface::setMagnet(const double b[3]) {
    gimic_set_magnet(b);
}

void GimicInterface::setSpin(char *s) {
    gimic_set_spin(s);
}

void GimicInterface::setScreening(double thrs) {
    gimic_set_screening(&thrs);
}

void GimicInterface::calcJTensor(const double r[3], double jt[9]) {
    gimic_calc_jtensor(r, jt);
}

void GimicInterface::calcJVector(const double r[3], double jv[3]) {
    gimic_calc_jvector(r, jv);
}

void GimicInterface::calcDivJ(const double r[3], double *dj) {
    gimic_calc_divj(r, dj);
}

void GimicInterface::calcModJ(const double r[3], double *mj) {
    gimic_calc_modj(r, mj);
}

void GimicInterface::calcAnapole(const double r[3], double *aj) {
    gimic_calc_anapole(r, aj);
}

