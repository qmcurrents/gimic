#ifndef MRENV_H
#define MRENV_H

class MREnv {
public:
    static void initializeMRCPP(int argc, char **argv, const char *fname = 0);
private:
    static void initializeTrees(int order, int max_depth,
                                double rel_prec, int polytype);
};
#endif // MRENV_H
