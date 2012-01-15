/**
 * \author Jonas Juselius <jonas.juselius@uit.no>, University of Tromsø
 */

#ifndef CURRENTFUNC_H_
#define CURRENTFUNC_H_

#include <iostream>
#include <cmath>

#include "BasicFunction.h"
#include "SeparableFunction.h"
#include "GimicFunc.h"


class CurrentFunc: public GimicFunc {
public:
        CurrentFunc(GimicInterface &g, int comp);
        virtual ~CurrentFunc() {}
        virtual double evalf(const double *r) const;
protected:
    int component;
};

#endif /* CURRENTFUNC_H_ */
