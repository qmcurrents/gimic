/**
 *
 */

#include "ModulusFunc.h"

using namespace std;

ModulusFunc::ModulusFunc(GimicInterface &g) : GimicFunc(g) {
}

double ModulusFunc::evalf(const double *r) const {
    double *j;
    gimic->calcModJ(r, j);
    return *j;
}
