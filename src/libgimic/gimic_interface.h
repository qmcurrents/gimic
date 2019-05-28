#ifndef GIMIC_CIF_H
#define GIMIC_CIF_H

#include "config.h"

#ifdef __cplusplus
extern "C" {
#endif
	void gimic_init(const char *, const char *);
    void gimic_finalize();
    void gimic_set_uhf(int *);
    void gimic_set_magnet(const double *);
    void gimic_set_spin(const char *);
    void gimic_set_screening(const double *);
    void gimic_calc_jtensor(const double *, double *);
    void gimic_calc_jvector(const double *, double *);
    void gimic_calc_modj(const double *, double *);
#ifdef __cplusplus
}
#endif

#endif
