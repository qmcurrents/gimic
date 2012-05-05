#ifndef GAUSSPOINTS_H
#define GAUSSPOINTS_H

#include "config.h"
#include "fcmangle.h"

#define gimic_get_gauss_points FC_MODULE_(\
        gausspoints,\
        gimic_get_gauss_points,\
        GAUSSPOINTS,\
        GIMIC_GET_GAUSS_POINTS)

void gimic_get_gauss_points(double *, double *, int *, int *, 
        double *, double *);
void mkgausspoints(double *, double *, int *, int *, double *, double *);

#endif
