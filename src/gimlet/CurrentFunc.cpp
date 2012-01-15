/**
 *
 */

#include "NodeIndex.h"
#include "ProjectedNode.h"
#include "CurrentFunc.h"

using namespace Eigen;
using namespace std;

CurrentFunc::CurrentFunc(GimicInterface &g, int comp) : GimicFunc(g) {
    assert(comp >= 0 & comp <= 2);
    component = comp;
}

double CurrentFunc::evalf(const double *r) const {
    Vector3d j;
    gimic->calcJVector(r, j.data());
    return j(component);
}
