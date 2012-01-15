/**
 *
 */

#include "NodeIndex.h"
#include "ProjectedNode.h"
#include "CurrentFunc.h"

using namespace Eigen;
using namespace std;

GimicFunc::GimicFunc(GimicInterface &g) {
    this->gimic = &g;
}

bool GimicFunc::isVisibleAtScale(int scale, int nQuadPts) const {
    int visibleScale = 1;
    if (scale < visibleScale) {
        return false;
    }
    return true;
}

