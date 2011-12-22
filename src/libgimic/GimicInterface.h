#ifndef GIMIC_INTERFACE_H
#define GIMIC_INTERFACE_H

class GimicInterface {
    public:
        GimicInterface(const char *mol, const char *xdens);
        virtual ~GimicInterface();
};
#endif
