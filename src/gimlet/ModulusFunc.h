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

        friend std::ostream& operator<<(std::ostream &o,
                                        const ModulusFunc &func)
	{
                o << "Not implemented yet." << std::endl;
		return o;
	}
protected:
};

#endif /* MODULUSFUNC_H_ */
