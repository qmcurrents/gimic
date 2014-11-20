#ifndef GAUSSPOINTS_H
#define GAUSSPOINTS_H

#include "config.h"

#ifdef __cplusplus
extern "C" {
#endif

void gimic_get_gauss_points(double *, double *, int *, int *, 
        double *, double *);
void mkgausspoints(double *, double *, int *, int *, double *, double *);

#ifdef __cplusplus
}
#endif

#endif
