#include "GimicInterface.h"
#include "gimic_cif.h"

GimicInterface::GimicInterface(const char *mol, const char *xdens) {
    gimic_init(mol, xdens);
}
GimicInterface::~GimicInterface() {
    gimic_finalize();
}
