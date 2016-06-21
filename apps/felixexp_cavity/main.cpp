#include <SOSkMC.h>

#include "../apputils.h"

int main(int argv, char** argc)
{
    rng.initialize(time(nullptr));

    string cfgName = getCfgName(argv, argc, "felixexp_cavity");

    Config cfg;
    cfg.readFile(cfgName.c_str());

    const Setting &cfgRoot = cfg.getRoot();

    const string &path = getSetting<string>(cfgRoot, "path") + "/";

    const uint &L = getSetting<uint>(cfgRoot, "L");
    const uint &W = getSetting<uint>(cfgRoot, "W");

    const double &height = getSetting<double>(cfgRoot, "height");
    const double &alpha = getSetting<double>(cfgRoot, "alpha");

    const double &flux = getSetting<double>(cfgRoot, "flux");

    const uint &boundaryDepth = getSetting<uint>(cfgRoot, "boundaryDepth");

    const uint &nCycles = getSetting<uint>(cfgRoot, "nCycles");
    const uint &interval = getSetting<uint>(cfgRoot, "interval");
    const uint &output = getSetting<uint>(cfgRoot, "output");

    SOSSolver solver(L, W, alpha, 0, true);

    const uint depositionBoxHalfSize = 3;
    const double maxDt = 0.01;

    RadialFirstPassage diff(solver, maxDt, depositionBoxHalfSize, getMFPTConstant(height, alpha, 0));
    FixedSurface confSurface(solver, height);

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

    AverageHeight h(solver);
    GrowthSpeed v(solver);

    if (output)
    {
        lattice.addEvent(h);
        lattice.addEvent(v);
        lattice.addEvent(dumper);
    }

    else
    {
        lattice.addEvent(storeHeights);
        lattice.addEvent(storeParticles);
    }

    initializeSurface(solver, "flat");

    const uint stepSize = 2*boundaryDepth;

    for (uint x = 0; x < L; ++x)
    {
        for (uint y = 0; y < stepSize; ++y)
        {
            solver.setHeight(x, y, 1, false);
        }
    }

    //app specifics
    const uint stepDim = 1;
    TrackAreaAverage y0AreaTracker(solver, 0, stepDim, boundaryDepth);

    Periodic x0(L, Boundary::orientations::FIRST);
    Periodic x1(L, Boundary::orientations::LAST);

    BlockByTracker y0(y0AreaTracker);
    Reflecting y1(W-1, Boundary::orientations::LAST);

    solver.setBoundaries({{&x0, &x1}, {&y0, &y1}});

    solver.registerPreNeighborObserver(&y0AreaTracker);

    PartialBoundaryNeighbors y0PartialNeighbors(solver, y0AreaTracker);
    NoBoundaryNeighbors y1NoNeighbors(solver, 0, y1.location(), stepDim);

    solver.addLocalPotential(&y0PartialNeighbors);
    solver.addLocalPotential(&y1NoNeighbors);

    solver.addFluxBoundary(stepDim, Boundary::orientations::FIRST, flux);

    lattice.eventLoop(nCycles);

    simRoot["alpha"] = alpha;
    simRoot["flux"] = flux;
    simRoot["height"] = height;
    simRoot["time"] = colvec(lattice.storedEventValues().col(0));

    return 0;
}
