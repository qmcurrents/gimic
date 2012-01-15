#ifndef GIMIC_CIF_H
#define GIMIC_CIF_H

#include "config.h"

#define gimic_init FC_FUNC_(gimic_init, GIMIC_INIT)
#define gimic_finalize FC_FUNC_(gimic_finalize, GIMIC_FINALIZE)
#define gimic_set_uhf FC_FUNC_(gimic_set_uhf, GIMIC_SET_UHF)
#define gimic_set_magnet FC_FUNC_(gimic_set_magnet, GIMIC_SET_MAGNET)
#define gimic_set_spin FC_FUNC_(gimic_set_spin, GIMIC_SET_SPIN)
#define gimic_set_screening FC_FUNC_(gimic_set_screening, GIMIC_SET_SCREENING)
#define gimic_calc_jtensor FC_FUNC_(gimic_calc_jtensor, GIMIC_CALC_JTENSOR)
#define gimic_calc_jvector FC_FUNC_(gimic_calc_jvector, GIMIC_CALC_JVECTOR)
#define gimic_calc_divj FC_FUNC_(gimic_calc_divj, GIMIC_CALC_DIVJ)
#define gimic_calc_edens FC_FUNC_(gimic_calc_edens, GIMIC_CALC_EDENS)

extern "C" {
	void gimic_init(const char *, const char *);
    void gimic_finalize();
    void gimic_set_uhf(int *); 
    void gimic_set_magnet(const double *); 
    void gimic_set_spin(char *); 
    void gimic_set_screening(double *); 
    void gimic_calc_jtensor(const double *, double *); 
    void gimic_calc_jvector(const double *, double *); 
    void gimic_calc_divj(const double *, double *); 
    void gimic_calc_edens(const double *, double *); 
}

#endif
