#include <SOSkMC.h>

#include "../apputils.h"

#include <chrono>

using ignis::Lattice;

int main(int argv, char** argc)
{
    string cfgName = getCfgName(argv, argc, "diff_model_times");

    Config cfg;
    cfg.readFile(cfgName.c_str());

    const Setting & cfgRoot = cfg.getRoot();

    const string path = getSetting<string>(cfgRoot, "path") + "/";

    const uint &L = getSetting<uint>(cfgRoot, "L");
    const uint &W = getSetting<uint>(cfgRoot, "W");

    const double &alpha = getSetting<double>(cfgRoot, "alpha");
    const double &height = getSetting<double>(cfgRoot, "height");

    const uint &diffusionType = getSetting<uint>(cfgRoot, "diffusionType");

    const uint &nCycles = getSetting<uint>(cfgRoot, "nCycles");
    const uint &nSurfaceEvents = getSetting<uint>(cfgRoot, "nSurfaceEvents");

    const uint &stdoutInterval = getSetting<uint>(cfgRoot, "stdoutInterval");

    rng.initialize((time(nullptr) % 1000000) + getProc(argv, argc));

    SOSSolver solver(L, W, alpha, 0, true);

    setBoundariesFromIDs(&solver, {0,0,0,0}, L, W);

    Diffusion* diff;

    switch (diffusionType) {
    case 1:
        diff = new ConfinedConstantConcentration(solver);
        break;
    case 2:
        diff = new LatticeDiffusion(solver);
        break;
    case 3:
        diff = new RadialFirstPassage(solver, 0.01, 3, getMFPTConstant(height, alpha, 0));
        break;
    case 4:
        diff = new AStarFirstPassage(solver, 0.01, 3, getMFPTConstant(height, alpha, 1));
        break;
    default:
        cout << "invalid diffusion type " << diffusionType << endl;
        return 1;
        break;
    }

    ConstantConfinement conf(solver, height);

    Lattice lattice;

    lattice.addEvent(solver);
    lattice.addEvent(diff);
    lattice.addEvent(conf);

    SurfaceCounter counter(nSurfaceEvents);
    solver.registerObserver(&counter);
    lattice.addEvent(counter);

    H5Wrapper::Root h5root(path +  addProcEnding(argv, argc, "diff_model_times", "h5"));
    H5Wrapper::Member &simRoot = setuph5(h5root, getProc(argv, argc), L, W);

    lattice.enableOutput(true, stdoutInterval);
    lattice.enableProgressReport();

    initializeSurface(solver, "random");

    /*
     *
     */

    auto t0 = chrono::steady_clock::now();
    lattice.eventLoop(nCycles);
    auto t1 = chrono::steady_clock::now();

    auto totalTime = t1 - t0;
    double cpuTime = chrono::duration<double, milli> (totalTime).count();

    /*
     *
     */

    simRoot["alpha"] = alpha;
    simRoot["height"] = height;
    simRoot["diffusionType"] = diffusionType;
    simRoot["cpuTime"] = cpuTime;

    return 0;
}
