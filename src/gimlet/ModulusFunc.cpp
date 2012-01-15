/**
 *
 */

#include "NodeIndex.h"
#include "ProjectedNode.h"
#include "ModulusFunc.h"

using namespace Eigen;
using namespace std;

ModulusFunc::ModulusFunc(GimicInterface &g) : GimicFunc(g) {
}

double ModulusFunc::evalf(const double *r) const {
    Vector3d j;
    gimic->calc_jvector(r, j.data());
    return j(0);
}
