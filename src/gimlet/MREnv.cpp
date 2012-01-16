#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>

#include "MREnv.h"
#include "Getkw.h"
#include "TelePrompter.h"
#include "LegendreBasis.h"
#include "InterpolatingBasis.h"
#include "QuadratureCache.h"
#include "FunctionTree.h"
#include "parallel.h"

using namespace std;

void MREnv::initializeMRCPP(int argc, char **argv, const char *fname) {
    omp_set_dynamic(0);
    mpi::communicator world;

    int nLocales = world.size();
    int nThreads = omp_get_max_threads();

    Eigen::internal::setNbThreads(1);
    const char *infile = 0;
    if (argc == 1) {
            infile = "STDIN";
    } else if (argc == 2) {
            infile = argv[1];
    } else {
            MSG_ERROR("Ivalid number of arguments!");
    }
    Getkw Input = Getkw(infile, false, true);

    int printLevel = Input.get<int>("printlevel");
    bool teletype = Input.get<bool>("teletype");

    if (fname != 0) {
            TelePrompter::init(printLevel, teletype, fname);
    } else {
            TelePrompter::init(printLevel, teletype, "GIMLET");
    }

    println(0, endl << endl);
    println(0, "*** Using the wavlet code from the MRCPP project by:");
    println(0, "***    Jonas Juselius   <jonas.juselius@uit.no> ");
    println(0, "***    Stig Rune Jensen <stig.r.jensen@uit.no>  ");
    println(0, "***    Luca Frediani    <luca.frediani@uit.no>  ");

    if (nLocales > 1 or nThreads > 1) {
        println(1, "+++ Parallel execution: ");
        println(1, "  Num MPI hosts  : " << nLocales);
        println(1, "  Threads/host   : " << nThreads);
        println(1, "  Total CPUs     : " << nLocales * nThreads);
        println(1, "");
    } else {
        println(1, "+++ Serial execution" << endl);
    }

    // initialize QuadratureCache globally to [0.1]
    getQuadratureCache(qCache);
    qCache.setBounds(0.0, 1.0);

    //Initialize world
    int order = Input.get<int>("Gimlet.order");
    int max_depth = Input.get<int>("Gimlet.max_depth");
    double rel_prec = Input.get<double>("Gimlet.rel_prec");
    string wlet = Input.get<string>("Gimlet.wavelet");

    int rootScale = Input.get<int>("Gimlet.initial_scale");
    const vector<int> &nbox = Input.getIntVec("Gimlet.boxes");
    const vector<int> &transl = Input.getIntVec("Gimlet.translation");
    const vector<double> &origin = Input.getDblVec("Gimlet.origin");

    BoundingBox<1>::setWorldBox(rootScale, transl.data(), nbox.data(), origin.data());
    BoundingBox<2>::setWorldBox(rootScale, transl.data(), nbox.data(), origin.data());
    BoundingBox<3>::setWorldBox(rootScale, transl.data(), nbox.data(), origin.data());

    const BoundingBox<3> &worldbox = BoundingBox<3>::getWorldBox();
    println(1, worldbox);

    int polytype;
    if (wlet == "I") {
        polytype = Interpol;
    } else {
        polytype = Legendre;
    }
    println(1, "*Default parameters:");
    println(1, "  Debug level  :     " << printLevel);
    println(1, "  Default order:     " << order);
    println(1, "  Default max depth: " << max_depth);
    println(1, "  Default precision: " << rel_prec);
    printout(1, "  Default polynomial type: ");
    if (polytype == Interpol) println(1, "Interpolating");
    if (polytype == Legendre) println(1, "Legendre");
    println(1, endl);
    println(1, endl);

    initializeTrees(order, max_depth, rel_prec, polytype);
}

void MREnv::initializeTrees(int order, int max_depth, double rel_prec,
                            int polytype) {
    FunctionTree<1>::setDefaultOrder(order);
    FunctionTree<1>::setDefaultMaxDepth(max_depth);
    FunctionTree<1>::setDefaultPrecision(rel_prec);
    FunctionTree<1>::setDefaultScalingType(polytype);

    FunctionTree<2>::setDefaultOrder(order);
    FunctionTree<2>::setDefaultMaxDepth(max_depth);
    FunctionTree<2>::setDefaultPrecision(rel_prec);
    FunctionTree<2>::setDefaultScalingType(polytype);

    FunctionTree<3>::setDefaultOrder(order);
    FunctionTree<3>::setDefaultMaxDepth(max_depth);
    FunctionTree<3>::setDefaultPrecision(rel_prec);
    FunctionTree<3>::setDefaultScalingType(polytype);

}

