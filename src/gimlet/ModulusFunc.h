/**
 * \author Jonas Juselius <jonas.juselius@uit.no>, University of Tromsø
 */

#ifndef MODULUSFUNC_H_
#define MODULUSFUNC_H_

#include <iostream>
#include <cmath>

#include "BasicFunction.h"
#include "SeparableFunction.h"
#include "GimicFunc.h"


class ModulusFunc: public GimicFunc {
public:
        ModulusFunc(GimicInterface &g);
        virtual ~ModulusFunc() {}
        virtual double evalf(const double *r) const;
};

#endif /* MODULUSFUNC_H_ */
