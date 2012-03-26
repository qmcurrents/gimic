#ifndef GIMIC_CIF_H
#define GIMIC_CIF_H

#include "config.h"
#include "fcmangle.h"

#define gimic_init FC_MODULE_(\
        gimic_interface,\
        gimic_init,\
        GIMIC_INTERFACE,\
        GIMIC_INIT)
#define gimic_finalize FC_MODULE_(\
        gimic_interface,\
        gimic_finalize,\
        GIMIC_INTERFACE,\
        GIMIC_FINALIZE)
#define gimic_set_uhf FC_MODULE_(\
        gimic_interface,\
        gimic_set_uhf,\
        GIMIC_INTERFACE,\
        GIMIC_SET_UHF)
#define gimic_set_magnet FC_MODULE_(\
        gimic_interface,\
        gimic_set_magnet,\
        GIMIC_INTERFACE,\
        GIMIC_SET_MAGNET)
#define gimic_set_spin FC_MODULE_(\
        gimic_interface,\
        gimic_set_spin,\
        GIMIC_INTERFACE,\
        GIMIC_SET_SPIN)
#define gimic_set_screening FC_MODULE_(\
        gimic_interface,\
        gimic_set_screening,\
        GIMIC_INTERFACE,\
        GIMIC_SET_SCREENING)
#define gimic_calc_jtensor FC_MODULE_(\
        gimic_interface,\
        gimic_calc_jtensor,\
        GIMIC_INTERFACE,\
        GIMIC_CALC_JTENSOR)
#define gimic_calc_jvector FC_MODULE_(\
        gimic_interface,\
        gimic_calc_jvector,\
        GIMIC_INTERFACE,\
        GIMIC_CALC_JVECTOR)
#define gimic_calc_divj FC_MODULE_(\
        gimic_interface,\
        gimic_calc_divj,\
        GIMIC_INTERFACE,\
        GIMIC_CALC_DIVJ)
#define gimic_calc_edens FC_MODULE_(\
        gimic_interface,\
        gimic_calc_edens,\
        GIMIC_INTERFACE,\
        GIMIC_CALC_EDENS)

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
