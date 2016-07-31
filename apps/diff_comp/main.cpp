#include <SOSkMC.h>

#include "../apputils.h"

int main(int argv, char** argc)
{
    rng.initialize(time(nullptr) + getProc(argv, argc));

    string cfgName = getCfgName(argv, argc, "diff_comp");

    Config cfg;
    cfg.readFile(cfgName.c_str());

    const Setting &cfgRoot = cfg.getRoot();

    const string &path = getSetting<string>(cfgRoot, "path") + "/";

    const uint &L = getSetting<uint>(cfgRoot, "L");
    const uint &W = getSetting<uint>(cfgRoot, "W");

    const double &height = getSetting<double>(cfgRoot, "height");
    const double &alpha = getSetting<double>(cfgRoot, "alpha");

    const double &flux = getSetting<double>(cfgRoot, "flux");

    const bool reflecting = getSetting<uint>(cfgRoot, "reflecting") == 1;

    const uint &nCycles = getSetting<uint>(cfgRoot, "nCycles");
    const uint &interval = getSetting<uint>(cfgRoot, "interval");
    const uint &output = getSetting<uint>(cfgRoot, "output");

    SOSSolver solver(L, W, alpha, 0);

    const uint depositionBoxHalfSize = 3;
    const double maxDt = 0.01;

    RadialFirstPassage diff(solver, maxDt, depositionBoxHalfSize, getMFPTConstant(height, alpha, 0));
    FixedSurface confSurface(solver, height);

    Time time(solver);

    DumpSystem dumper(solver, interval, path, getTail(argv, argc));

    H5Wrapper::Root h5root(path +  addProcEnding(argv, argc, "diff_comp", "h5"));
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

    lattice.addEvent(storeHeights);
    lattice.addEvent(storeParticles);

    initializeSurface(solver, "thermalized");

    Periodic x0(L, Boundary::orientations::FIRST);
    Periodic x1(L, Boundary::orientations::LAST);

    ConstantHeight y0(0, 0, Boundary::orientations::FIRST);
    solver.addFluxBoundary(1, Boundary::orientations::FIRST, flux);

    Boundary *y1;
    NoBoundaryNeighbors y1NoNeighbors(solver, 0, W-1, 1);
    if (reflecting)
    {
        y1 = new Reflecting(W-1, Boundary::orientations::LAST);
        solver.addLocalPotential(&y1NoNeighbors);
    }

    else
    {
        y1 = new ConstantHeight(0, W-1, Boundary::orientations::LAST);
        solver.addFluxBoundary(1, Boundary::orientations::LAST, 1.0);
    }

    solver.setBoundaries({{&x0, &x1}, {&y0, y1}});

    lattice.eventLoop(nCycles);

    simRoot["alpha"] = alpha;
    simRoot["flux"] = flux;
    simRoot["reflecting"] = reflecting;
    simRoot["height"] = height;
    simRoot["time"] = colvec(lattice.storedEventValues().col(0));


    delete y1;

    return 0;
}
