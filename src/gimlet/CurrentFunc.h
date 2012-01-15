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
        CurrentFunc(GimicInterface &g);
        virtual ~CurrentFunc() {}

        virtual double evalf(const double *r) const;

        friend std::ostream& operator<<(std::ostream &o,
                                        const CurrentFunc &func)
	{
                o << "Not implemented yet." << std::endl;
		return o;
	}
protected:

};

#endif /* CURRENTFUNC_H_ */
