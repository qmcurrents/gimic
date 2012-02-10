#ifndef GIMIC_CIF_H
#define GIMIC_CIF_H

#include "config.h"
#include "fcmangle.h"

#define gimic_init FC_GLOBAL(gimic_init, GIMIC_INIT)
#define gimic_finalize FC_GLOBAL(gimic_finalize, GIMIC_FINALIZE)
#define gimic_set_uhf FC_GLOBAL(gimic_set_uhf, GIMIC_SET_UHF)
#define gimic_set_magnet FC_GLOBAL(gimic_set_magnet, GIMIC_SET_MAGNET)
#define gimic_set_spin FC_GLOBAL(gimic_set_spin, GIMIC_SET_SPIN)
#define gimic_set_screening FC_GLOBAL(gimic_set_screening, GIMIC_SET_SCREENING)
#define gimic_calc_jtensor FC_GLOBAL(gimic_calc_jtensor, GIMIC_CALC_JTENSOR)
#define gimic_calc_jvector FC_GLOBAL(gimic_calc_jvector, GIMIC_CALC_JVECTOR)
#define gimic_calc_divj FC_GLOBAL(gimic_calc_divj, GIMIC_CALC_DIVJ)
#define gimic_calc_edens FC_GLOBAL(gimic_calc_edens, GIMIC_CALC_EDENS)

extern "C" {
	void gimic_init(const char *, const char *);
    void gimic_finalize();
    void gimic_set_uhf(int *); 
    void gimic_set_magnet(double *); 
    void gimic_set_spin(char *); 
    void gimic_set_screening(double *); 
    void gimic_calc_jtensor(double *, double *); 
    void gimic_calc_jvector(double *, double *); 
    void gimic_calc_divj(double *, double *); 
    void gimic_calc_edens(double *, double *); 
}

#endif
