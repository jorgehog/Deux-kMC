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

    ReflConstantHybrid topBoundary(solver, 0, Boundary::orientations::LAST, 1);
    ReflConstantHybrid rightBoundary(solver, 0, Boundary::orientations::LAST, 0);

    AverageHeightLineBoundaryRefl leftBoundary(solver, Boundary::orientations::FIRST, 0, boundaryDepth);
    AverageHeightLineBoundaryOpen bottomBoundary(solver, Boundary::orientations::FIRST, 1, boundaryDepth);

    solver.setBoundaries({{&leftBoundary, &rightBoundary}, {&bottomBoundary, &topBoundary}});

    ConcentrationProfile concProfile(solver, [&W, &cBath] (const uint x, const uint y)
    {
        (void) x;
        return cBath - (cBath - 1)*y/(W-1);
    });

    ConstantConfinement confinement(solver, h0);

    Lattice lattice;
    lattice.enableProgressReport();

    lattice.addEvent(solver);
    lattice.addEvent(confinement);
    lattice.addEvent(concProfile);

    DumpSystem dumper(solver, interval, path);

    lattice.addEvent(dumper);

        for (uint x = 0; x < size; ++x)
        {
            const uint ym = sqrt(size*size - x*x);

            for (uint y = 0; y < ym; ++y)
            {
                solver.setHeight(x, y, 1, false);
            }
        }

    //
    lattice.eventLoop(nCycles);
    //

    H5Wrapper::Root h5root(path +  addProcEnding(argv, argc, "felixexp", "h5"));
    H5Wrapper::Member &simRoot = setuph5(h5root, getProc(argv, argc), L, W);

    simRoot["L"] = L;
    simRoot["W"] = W;
    simRoot["alpha"] = alpha;
    simRoot["h0"] = h0;

    simRoot["cBath"] = cBath;

    simRoot[""
            ""
            ""
            ""
            "e"] = halfSize;
    simRoot["maxdt"] = maxdt;

    return 0;
}
