/**
 * \author Jonas Juselius <jonas.juselius@uit.no>, University of Tromsø
 */

#ifndef GIMICFUNC_H_
#define GIMICTFUNC_H_

#include <iostream>
#include <cmath>

#include "BasicFunction.h"
#include "RepresentableFunction.h"
#include "GimicInterface.h"


class GimicFunc: public RepresentableFunction<3> {
public:
        GimicFunc(GimicInterface &g) {
             this->gimic = &g;
        }
        virtual ~GimicFunc() {}
        virtual double evalf(const double *r) const = 0;
protected:
        GimicInterface *gimic;

        bool isVisibleAtScale(int scale, int nQuadPts) const {
            int visibleScale = 1;
            if (scale < visibleScale) {
                return false;
            }
            return true;
        }
};

#endif /* GIMICFUNC_H_ */
