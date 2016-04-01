#include <SOSkMC.h>

#include "../apputils.h"

#include "extraneighbor.h"
#include "rdlextraneighborsurface.h"

using ignis::Lattice;

class Coverage : public SOSEvent
{
public:
    Coverage(const SOSSolver &solver) :
        SOSEvent(solver, "Coverage", "", true, true)
    {

    }

    // Event interface
public:
    void execute()
    {
        uint coverage = 0;
        const double &hl = solver().confiningSurfaceEvent().height();

        for (uint x = 0; x < solver().length(); ++x)
        {
            for (uint y = 0; y < solver().width(); ++y)
            {
                if (hl - solver().height(x, y) < 2)
                {
                    coverage++;
                }
            }
        }

        setValue(coverage);
    }
};

int main(int argv, char** argc)
{
    string cfgName = getCfgName(argv, argc, "extraneighbor");

    Config cfg;
    cfg.readFile(cfgName.c_str());

    const Setting & cfgRoot = cfg.getRoot();

    const string path = getSetting<string>(cfgRoot, "path") + "/";
    const double Pl = getSetting<double>(cfgRoot, "Pl");
    const double alpha = getSetting<double>(cfgRoot, "alpha");
    const double ld = getSetting<double>(cfgRoot, "ld");
    const double s0 = getSetting<double>(cfgRoot, "s0");
    const uint nCycles = getSetting<uint>(cfgRoot, "nCycles");
    const uint interval = getSetting<uint>(cfgRoot, "interval");
    const uint output = getSetting<uint>(cfgRoot, "output");

    rng.initialize(time(nullptr));
//    rng.initialize(10000);

    const uint L = 30;
    const uint W = 30;

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

    Lattice lattice;

    lattice.addEvent(solver);
    lattice.addEvent(rdlSurface);
    lattice.addEvent(constantConcentration);

    Coverage coverage(solver);
    lattice.addEvent(coverage);

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
                                    output,
                                    "ignisSOS.ign",
                                    "/tmp",
                                    interval);

    initializeSurface(solver, "random");
    rdlSurface.setHeight(rdlSurface.getRdlEquilibrium());

    lattice.eventLoop(nCycles);

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
    simRoot["coverage"] = colvec(lattice.storedEventValues().col(0));

    return 0;
}
