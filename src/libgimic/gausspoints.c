#include <stdio.h>
#include "gausspoints.h"

void mkgausspoints(double *a, double *b, int *npts, double *pts, double *wgts) {
    gimic_get_gauss_points(a, b, npts, pts, wgts);
}

