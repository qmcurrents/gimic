#include <iostream>
#include <Eigen/Core>
#include "GimicInterface.h"

using namespace std;
using namespace Eigen;

int main() {
    cout << "GIMLET starting." << endl;
    Vector3d b, r, j;
    b << 0.0, 0.0 , 1.0;
    r << 0.0, 0.1, 0.1;
    j << 0.0, 0.0, 0.0;

    GimicInterface gimic("mol", "xdens");
    gimic.set_magnet(b.data());
//    gimic.set_spin("total");
    gimic.set_uhf(0);
    gimic.calc_jvector(r.data(), j.data());
    cout << "Magnet: " << endl << b << endl;
    cout << "Coord: " << endl << r << endl;
    cout << "Current: " << endl << j << endl;

    return 0;
}
