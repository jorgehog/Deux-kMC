#include <SOSkMC.h>

#include "../apputils.h"

using ignis::Lattice;

int main(int argv, char** argc)
{
    string cfgName = getCfgName(argv, argc, "felixexp_cavity");

    Config cfg;
    cfg.readFile(cfgName.c_str());

    const Setting & cfgRoot = cfg.getRoot();

    const string &path = getSetting<string>(cfgRoot, "path") + "/";

    const uint &L = getSetting<uint>(cfgRoot, "L");
    const uint &W = getSetting<uint>(cfgRoot, "W");

    const double &alpha = getSetting<double>(cfgRoot, "alpha");
    const double &omega = getSetting<double>(cfgRoot, "omega");
    const double &height = getSetting<double>(cfgRoot, "height");

    const uint &boundaryDepth = getSetting<uint>(cfgRoot, "boundaryDepth");

    const uint &nCycles = getSetting<uint>(cfgRoot, "nCycles");
    const uint &interval = getSetting<uint>(cfgRoot, "interval");
    const uint &output = getSetting<uint>(cfgRoot, "output");

    rng.initialize(time(nullptr));
    const string initType = "flat";

    const double gamma = log(1 + omega);

    SOSSolver solver(L, W, alpha, gamma, true);

    ReflectingSurfaceOpenSolution x0(0, Boundary::orientations::FIRST);
    Reflecting x1(L-1, Boundary::orientations::LAST);

    Periodic y0(W, Boundary::orientations::FIRST);
    Periodic y1(W, Boundary::orientations::LAST);

    solver.setBoundaries({{&x0, &x1}, {&y0, &y1}});

    TrackAreaAverage inletAreaTracker(solver, 0, 0, boundaryDepth, true);
    TrackAreaAverage cavityAreaTracker(solver, L-1, 0, boundaryDepth, true);

    solver.registerObserver(&inletAreaTracker);
    solver.registerObserver(&cavityAreaTracker);

    PartialBoundaryNeighbors inletPartialNeighbors(solver, inletAreaTracker);
    PartialBoundaryNeighbors cavityPartialNeighbors(solver, cavityAreaTracker);

    solver.addLocalPotential(&inletPartialNeighbors);
    solver.addLocalPotential(&cavityPartialNeighbors);

    const uint depositionBoxHalfSize = 3;
    const double maxDt = 0.01;

    RadialFirstPassage diff(solver, maxDt, depositionBoxHalfSize, getMFPTConstant(height, alpha, 0));
    FixedSurface confSurface(solver, height);

    solver.addConcentrationBoundary(0, Boundary::orientations::FIRST, omega);

    Time time(solver);

    DumpSystem dumper(solver, interval, path, getTail(argv, argc));

    H5Wrapper::Root h5root(path +  addProcEnding(argv, argc, "felix_cavity", "h5"));
    H5Wrapper::Member &simRoot = setuph5(h5root, getProc(argv, argc), L, W);

    StoreHeights storeHeights(solver, interval, simRoot);
    StoreParticles storeParticles(solver, diff, interval, simRoot);

    Lattice lattice;

    lattice.enableOutput(true, interval);
    lattice.enableProgressReport();
    lattice.enableEventValueStorage(true,
                                    output,
                                    "ignisSOS.ign",
                                    "/tmp",
                                    interval);

    lattice.addEvent(solver);
    lattice.addEvent(confSurface);
    lattice.addEvent(diff);
    lattice.addEvent(time);

    if (output)
    {
        lattice.addEvent(dumper);
    }
    else
    {
        lattice.addEvent(storeHeights);
        lattice.addEvent(storeParticles);
    }

    initializeSurface(solver, initType, 3, 10000);

    lattice.eventLoop(nCycles);

    simRoot["alpha"] = alpha;
    simRoot["omega"] = omega;
    simRoot["height"] = height;
    simRoot["time"] = colvec(lattice.storedEventValues().col(0));

    return 0;
}
