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
        GimicFunc(GimicInterface &g);
        virtual ~GimicFunc() {}

        virtual double evalf(const double *r) const = 0;

        friend std::ostream& operator<<(std::ostream &o,
                                        const GimicFunc &func)
	{
                o << "Not implemented yet." << std::endl;
		return o;
	}
protected:
        GimicInterface *gimic;

        bool isVisibleAtScale(int scale, int nQuadPts) const;
};

#endif /* GIMICFUNC_H_ */
