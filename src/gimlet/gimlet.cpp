#include <iostream>
#include <Eigen/Core>
#include "GimicInterface.h"
#include "FunctionTree.h"
#include "ModulusFunc.h"
#include "parallel.h"
#include "TelePrompter.h"
#include "MREnv.h"

using namespace std;
using namespace Eigen;

int main(int argc, char **argv) {
    cout << "GIMLET starting." << endl;

    mpi::environment env(argc, argv);
    TelePrompter::init();

    MREnv::initializeMRCPP(argc, argv);

    Vector3d b, r, j;
    b << 0.0, 0.0 , 1.0;
    r << 0.0, 0.1, 0.1;
    j << 0.0, 0.0, 0.0;

    GimicInterface gimic("mol", "xdens");
    gimic.setMagnet(b.data());
//    gimic.setSpin("total");
    gimic.setUhf(0);
    gimic.calcJVector(r.data(), j.data());
    cout << "Magnet: " << endl << b << endl;
    cout << "Coord: " << endl << r << endl;
    cout << "Current: " << endl << j << endl;

    ModulusFunc modj(gimic);
    FunctionTree<3> ftree;
    ftree.projectFunction(modj);
    int integral = ftree.integrate();

    cout << ftree << endl;
    cout << "Integral = " << integral << endl;


    return 0;
}
