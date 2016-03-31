#include <SOSkMC.h>

#include "../apputils.h"

#include "extraneighbor.h"
#include "rdlextraneighborsurface.h"

using ignis::Lattice;

int main(int argv, char** argc)
{
    string cfgName = getCfgName(argv, argc, "extraneighbor");

    Config cfg;
    cfg.readFile(cfgName.c_str());

    const Setting & cfgRoot = cfg.getRoot();

    const string path = getSetting<string>(cfgRoot, "path") + "/";
    const double Pl = getSetting<double>(cfgRoot, "Pl");
    const double alpha = getSetting<double>(cfgRoot, "alpha");

    rng.initialize(time(nullptr));
//    rng.initialize(10000);

    const uint L = 50;
    const uint W = 50;

    const double ld = 5.0;
    const double s0 = 1.0;
    const double eTerm = 1-exp(-1/ld);

    const double gamma = alpha*Pl/(ld*eTerm);

    SOSSolver solver(L, W, alpha, gamma, true);

    setBoundariesFromIDs(&solver, {0,0,0,0}, L, W);

    RDLPotential rdlpotential(solver, s0, ld);
    solver.addLocalPotential(&rdlpotential);
    solver.registerObserver(&rdlpotential);

    ExtraNeighbor extraNeighbor(solver, s0/eTerm);
    solver.addLocalPotential(&extraNeighbor);


    RDLExtraNeighborSurface rdlSurface(solver, rdlpotential, extraNeighbor, Pl);
//    FixedSurface rdlSurface(solver, 10);

//    RadialFirstPassage constantConcentration(solver, 0.01, 3, getMFPTConstant(10, alpha, 0));
    ConfinedConstantConcentration constantConcentration(solver);

    const uint interval = 100;

    Lattice lattice;

    lattice.addEvent(solver);
    lattice.addEvent(rdlSurface);
    lattice.addEvent(constantConcentration);

//    DumpSystem dumper(solver, interval, path);
//    lattice.addEvent(dumper);

    AverageHeight avgH(solver);
    lattice.addEvent(avgH);

#ifndef NDEBUG
    RateChecker checker(solver);
    lattice.addEvent(checker);
#endif

    lattice.enableOutput(true, interval);
    lattice.enableProgressReport();
    lattice.enableEventValueStorage(true,
                                    true,
                                    "ignisSOS.ign",
                                    "/tmp",
                                    interval);

    initializeSurface(solver, "random");
    rdlSurface.setHeight(rdlSurface.getRdlEquilibrium());

    lattice.eventLoop(1000000);


    H5Wrapper::Root h5root(path +  addProcEnding(argv, argc, "extraneighbor", "h5"));

    stringstream sizeDesc;
    sizeDesc << L << "x" << W;
    H5Wrapper::Member &sizeRoot = h5root.addMember(sizeDesc.str());

    timeval tv;
    gettimeofday(&tv, nullptr);

    __int64_t run_ID = 1000*tv.tv_sec + tv.tv_usec/1000 + 10000000000000u*getProc(argv, argc);

    H5Wrapper::Member &simRoot = sizeRoot.addMember(run_ID);

    simRoot["alpha"] = alpha;
    simRoot["Pl"] = Pl;
    simRoot["h"] = rdlSurface.height();
    simRoot["heights"] = solver.heights();

    return 0;
}
