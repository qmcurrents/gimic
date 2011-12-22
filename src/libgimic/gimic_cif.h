#ifndef GIMIC_CIF_H
#define GIMIC_CIF_H

#include "config.h"

#define gimic_init FC_FUNC_(gimic_init, GIMIC_INIT)
#define gimic_finalize FC_FUNC_(gimic_finalize, GIMIC_FINALIZE)

extern "C" {
	void gimic_init(const char *);
	void gimic_finalize();
}

#endif
