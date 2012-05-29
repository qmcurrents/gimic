#ifndef GIMIC_CIF_H
#define GIMIC_CIF_H

#include "config.h"

#ifdef __cplusplus
extern "C" {
#endif
	void gimic_init(const char *, const char *);
    void gimic_finalize();
    void gimic_set_uhf(int *); 
    void gimic_set_magnet(double *); 
    void gimic_set_spin(const char *); 
    void gimic_set_screening(double *); 
    void gimic_calc_jtensor(double *, double *); 
    void gimic_calc_jvector(double *, double *); 
    void gimic_calc_divj(double *, double *); 
    void gimic_calc_edens(double *, double *); 
#ifdef __cplusplus
}
#endif

#endif
