#include <SOSkMC.h>

#include "../apputils.h"

int main(int argv, char **argc)
{
    rng.initialize(time(nullptr) + getProc(argv, argc));

    string cfgName = getCfgName(argv, argc, "felixexp");

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

    DumpSystem dumper(solver, interval, path);

    H5Wrapper::Root h5root(path +  addProcEnding(argv, argc, "felixexp", "h5"));
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

    if (output == 1)
    {
        lattice.addEvent(dumper);
    }

    else
    {
        lattice.addEvent(storeHeights);
        lattice.addEvent(storeParticles);
    }

    initializeSurface(solver, "flat");

    //app specifics

    const uint &size = getSetting<uint>(cfgRoot, "size");
    for (uint y = 0; y < size; ++y)
    {
        const uint xm = L/W*sqrt(size*size - y*y);

        for (uint x = 0; x < xm; ++x)
        {
            solver.setHeight(x, y, 1, false);
        }
    }

    Reflecting x0(0, Boundary::orientations::FIRST);
    Reflecting x1(L-1, Boundary::orientations::LAST);

    const uint stepDim = 1;
    TrackLineAverage x0LineTracker(solver, 0, 0, boundaryDepth);
    TrackLineAverage y0LineTracker(solver, 0, 1, boundaryDepth);

    BlockByTracker y0(y0LineTracker);
    ConstantHeight y1(0, W-1, Boundary::orientations::LAST);

    solver.setBoundaries({{&x0, &x1}, {&y0, &y1}});

    solver.registerPreNeighborObserver(&x0LineTracker);
    solver.registerPreNeighborObserver(&y0LineTracker);

    PartialBoundaryNeighbors x0PartialNeighbors(solver, x0LineTracker);
    NoBoundaryNeighbors x1NoNeighbors(solver, 0, x1.location(), 0);

    PartialBoundaryNeighbors y0PartialNeighbors(solver, y0LineTracker);

    solver.addLocalPotential(&x0PartialNeighbors);
    solver.addLocalPotential(&y0PartialNeighbors);
    solver.addLocalPotential(&x1NoNeighbors);

    const double eqFlux = 0.75;
    solver.addFluxBoundary(stepDim, Boundary::orientations::FIRST, flux);
    solver.addFluxBoundary(stepDim, Boundary::orientations::LAST, eqFlux);

    lattice.eventLoop(nCycles);

    simRoot["alpha"] = alpha;
    simRoot["flux"] = flux;
    simRoot["height"] = height;
    simRoot["time"] = colvec(lattice.storedEventValues().col(0));

    return 0;
}
