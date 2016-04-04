#include <SOSkMC.h>

#include "miscevents.h"

#include "../apputils.h"

#include <HDF5Wrapper/include/hdf5wrapper.h>
#include <libconfig_utils/libconfig_utils.h>

int main(int argv, char **argc)
{
    rng.initialize(time(nullptr));
//    rng.initialize(139444);

    string cfgName = getCfgName(argv, argc, "felixexp");

    Config cfg;
    cfg.readFile(cfgName.c_str());

    const Setting & cfgRoot = cfg.getRoot();

    const string path = getSetting<string>(cfgRoot, "path") + "/";

    const int halfSize = getSetting<int>(cfgRoot, "halfSize");
    const double maxdt = getSetting<double>(cfgRoot, "maxdt");

    const uint L = getSetting<uint>(cfgRoot, "L");
    const uint W = getSetting<uint>(cfgRoot, "W");
    const double h0 = getSetting<double>(cfgRoot, "h0");
    const double alpha = getSetting<double>(cfgRoot, "alpha");

    const double cBath = getSetting<double>(cfgRoot, "cBath");

    const uint size = getSetting<uint>(cfgRoot, "size");
    const uint interval = getSetting<uint>(cfgRoot, "interval");

    const uint boundaryDepth = getSetting<uint>(cfgRoot, "boundaryDepth");

    const uint nCycles = getSetting<uint>(cfgRoot, "nCycles");


    SOSSolver solver(L, W, alpha, 0, true);

    AverageHeightBoundary leftBoundary(solver, boundaryDepth, 0, L, W, Boundary::orientations::FIRST, 0);
    ConstantHeight rightBoundary(0, L-1, Boundary::orientations::LAST);

//    Periodic leftBoundary(L, Boundary::orientations::FIRST);
//    Periodic rightBoundary(L, Boundary::orientations::LAST);


//    LongestStripBoundary bottomBoundary(solver, 1, Boundary::orientations::FIRST);
//    LongestStripBoundary topBoundary(solver, 1, Boundary::orientations::LAST);
//    ConstantHeight bottomBoundary(0, 0, Boundary::orientations::FIRST);
    ConstantHeight topBoundary(0, W-1, Boundary::orientations::LAST);
    AverageHeightLineBoundary bottomBoundary(solver, Boundary::orientations::FIRST, 1, boundaryDepth);
//    AverageHeightLineBoundary topBoundary(solver, Boundary::orientations::LAST, 1, boundaryDepth);


    solver.setBoundaries({{&leftBoundary, &rightBoundary}, {&bottomBoundary, &topBoundary}});

    ConcentrationProfile concProfile(solver, [&W, &cBath] (const uint x, const uint y)
    {
        (void) x;
        return cBath - (cBath - 1)*y/(W-1);
    });

    ConstantConfinement confinement(solver, h0);

    InsertStep insertStep(solver, interval, size);

    Lattice lattice;
    lattice.enableProgressReport();

    lattice.addEvent(solver);
    //    lattice.addEvent(insertStep);
    lattice.addEvent(confinement);
    lattice.addEvent(concProfile);


    DumpSystem dumper(solver, 1000, path);

    lattice.addEvent(dumper);

        for (uint x = 0; x < size; ++x)
        {
//            const uint ym = sqrt(size*size - x*x);

//            for (uint y = 0; y < ym; ++y)
//            {
//                solver.setHeight(x, y, 1, false);
//                solver.setHeight(L-x-1, y, 1, false);
//            }
//            solver.freezeSurfaceParticle(x, 0);
//            solver.freezeSurfaceParticle(L-x-1, 0);

            for (uint y = 0; y < W; ++y)
            {
                solver.setHeight(x, y, 1, false);
            }
        }

    //
    lattice.eventLoop(nCycles);
    //


    H5Wrapper::Root h5root(path +  addProcEnding(argv, argc, "felixexp", "h5"));

    stringstream sizeDesc;
    sizeDesc << L << "x" << W;
    H5Wrapper::Member &sizeRoot = h5root.addMember(sizeDesc.str());

    timeval tv;
    gettimeofday(&tv, nullptr);

    __int64_t run_ID = 1000*tv.tv_sec + tv.tv_usec/1000 + 10000000000000u*getProc(argv, argc);

    H5Wrapper::Member &simRoot = sizeRoot.addMember(run_ID);

    simRoot["L"] = L;
    simRoot["W"] = W;
    simRoot["alpha"] = alpha;
    simRoot["h0"] = h0;

    simRoot["cBath"] = cBath;

    simRoot["halfSize"] = halfSize;
    simRoot["maxdt"] = maxdt;

    return 0;
}
