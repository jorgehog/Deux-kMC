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

    const double omegaEq = getSetting<double>(cfgRoot, "omegaEq");
    const double omega = getSetting<double>(cfgRoot, "omega");

    const uint size = getSetting<uint>(cfgRoot, "size");
    const uint boundaryDepth = getSetting<uint>(cfgRoot, "boundaryDepth");

    const uint nCycles = getSetting<uint>(cfgRoot, "nCycles");
    const uint interval = getSetting<uint>(cfgRoot, "interval");
    const uint output = getSetting<uint>(cfgRoot, "output");

    SOSSolver solver(L, W, alpha, 0, true);

    //    ReflConstantHybrid topBoundary(solver, 0, Boundary::orientations::LAST, 1);
    ConstantHeight topBoundary(0, W-1, Boundary::orientations::LAST);
    ReflConstantHybrid rightBoundary(solver, 0, Boundary::orientations::LAST, 0);

    AverageHeightLineBoundaryRefl leftBoundary(solver, Boundary::orientations::FIRST, 0, boundaryDepth);
    AverageHeightLineBoundaryOpen bottomBoundary(solver, Boundary::orientations::FIRST, 1, boundaryDepth);

    solver.setBoundaries({{&leftBoundary, &rightBoundary}, {&bottomBoundary, &topBoundary}});


    //    ConcentrationProfile concProfile(solver, [&W, &cBath] (const uint x, const uint y)
    //    {
    //        (void) x;
    //        return cBath - (cBath - 1)*y/(W-1);
    //    });

    RadialFirstPassage diff(solver, maxdt, halfSize, getMFPTConstant(h0, alpha, 0));

    ConstantConfinement confinement(solver, h0);

    solver.addConcentrationBoundary(1, Boundary::orientations::FIRST, omega);
    solver.addConcentrationBoundary(1, Boundary::orientations::LAST, omegaEq);

    Lattice lattice;
    lattice.enableProgressReport();

    lattice.addEvent(solver);
    lattice.addEvent(confinement);
    lattice.addEvent(diff);

    DumpSystem dumper(solver, interval, path);

    H5Wrapper::Root h5root(path +  addProcEnding(argv, argc, "felixexp", "h5"));
    H5Wrapper::Member &simRoot = setuph5(h5root, getProc(argv, argc), L, W);

    StoreHeights storeHeights(solver, interval, simRoot);
    StoreParticles storeParticles(solver, diff, interval, simRoot);

    if (output == 1)
    {
        lattice.addEvent(dumper);
    }

    else
    {
        lattice.addEvent(storeHeights);
        lattice.addEvent(storeParticles);
    }

    Time time(solver);
    lattice.addEvent(time);

    lattice.enableEventValueStorage(true,
                                    false,
                                    "",
                                    "",
                                    interval);

    for (uint y = 0; y < size; ++y)
    {
        const uint xm = L/W*sqrt(size*size - y*y);

        for (uint x = 0; x < xm; ++x)
        {
            solver.setHeight(x, y, 1, false);
        }
    }

    //
    lattice.eventLoop(nCycles);
    //

    simRoot["alpha"] = alpha;
    simRoot["h0"] = h0;

    simRoot["omega"] = omega;
    simRoot["omegaEq"] = omegaEq;

    simRoot["time"] = colvec(lattice.storedEventValues().col(0));

    return 0;
}
