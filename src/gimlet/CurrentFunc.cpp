/**
 *
 */

#include "NodeIndex.h"
#include "ProjectedNode.h"
#include "CurrentFunc.h"

using namespace Eigen;
using namespace std;

CurrentFunc::CurrentFunc(GimicInterface &g) : GimicFunc(g) {
}

double CurrentFunc::evalf(const double *r) const {
    Vector3d j;
    gimic->calc_jvector(r, j.data());
    return j(0);
}
