#include <SOSkMC.h>

#include "../apputils.h"

#include "extraneighbor.h"
#include "rdlextraneighborsurface.h"

using ignis::Lattice;

class DumpForceProfile : public SOSEvent
{
public:
    DumpForceProfile(const RDLExtraNeighborSurface &surface,
                     const uint interval) :
        SOSEvent(surface.solver(),
                 "DumpForceProfile"),
        m_surface(surface),
        m_interval(interval)
    {

    }

private:
    const RDLExtraNeighborSurface &m_surface;
    const uint m_interval;

    // Event interface
public:
    void execute()
    {
        if (cycle() % m_interval == 0)
        {
            m_surface.dumpProfile(cycle()/m_interval);
        }
    }
};

class Coverage : public SOSEvent
{
public:
    Coverage(const SOSSolver &solver,
             const bool dumpCoverage = false,
             const uint interval = 1,
             const uint nCycles = 0) :
        SOSEvent(solver, "Coverage", "", true, true),
        m_dumpCoverage(dumpCoverage),
        m_interval(interval),
        m_coverage(dumpCoverage ? solver.length() : 0,
                   dumpCoverage ? solver.width()  : 0,
                   dumpCoverage ? nCycles/interval : 0,
                   fill::zeros)

    {

    }

    const icube &coverage() const
    {
        return m_coverage;
    }

private:
    const bool m_dumpCoverage;
    const uint m_interval;

    icube m_coverage;

    // Event interface
public:
    void execute()
    {
        const bool isDumpCycle = m_dumpCoverage && (cycle() % m_interval == 0);
        const uint n = cycle()/m_interval;

        const double &hl = solver().confiningSurfaceEvent().height();

        uint coverage = 0;
        for (uint x = 0; x < solver().length(); ++x)
        {
            for (uint y = 0; y < solver().width(); ++y)
            {
                if (hl - solver().height(x, y) < 2)
                {
                    if (isDumpCycle)
                    {
                        m_coverage(x, y, n) = 1;
                    }

                    coverage++;
                }
            }
        }

        setValue(coverage);
    }
};

class DetectZeroCoverage : public SOSEvent
{
public:
    DetectZeroCoverage(const Coverage *coverage, const uint nBeforeBreak) :
        SOSEvent(coverage->solver(), "DetectZeroCoverage"),
        m_coverage(coverage),
        m_nBeforeBreak(nBeforeBreak)
    {
        setDependency(coverage);
    }

private:

    const Coverage *m_coverage;
    const uint m_nBeforeBreak;

    uint n;


    // Event interface
public:
    void execute()
    {
        if (n >= m_nBeforeBreak)
        {
            terminateLoop("The coverage is zero. Ending...");
        }

        if (m_coverage->value() == 0)
        {
            n++;
        }

        else
        {
            n = 0;
        }

    }
    void initialize()
    {
        n = 0;
    }
};

int main(int argv, char** argc)
{
    string cfgName = getCfgName(argv, argc, "extraneighbor");

    Config cfg;
    cfg.readFile(cfgName.c_str());

    const Setting & cfgRoot = cfg.getRoot();

    const string path = getSetting<string>(cfgRoot, "path") + "/";

    const uint &L = getSetting<uint>(cfgRoot, "L");
    const uint &W = getSetting<uint>(cfgRoot, "W");

    const double &F0 = getSetting<double>(cfgRoot, "F0");
    const double &alpha = getSetting<double>(cfgRoot, "alpha");
    const int &omegaSign = getSetting<int>(cfgRoot, "omegaSign");
    const double &omegaVal = getSetting<double>(cfgRoot, "omegaVal");

    const double &ld = getSetting<double>(cfgRoot, "ld");
    const double &s0 = getSetting<double>(cfgRoot, "s0");

    const uint &nCycles = getSetting<uint>(cfgRoot, "nCycles");
    const uint &nCyclesOmega = getSetting<uint>(cfgRoot, "nCyclesOmega");
    const uint &ignisOutput = getSetting<uint>(cfgRoot, "ignisOutput");
    const uint &stdoutInterval = getSetting<uint>(cfgRoot, "stdoutInterval");
    const uint &interval = getSetting<uint>(cfgRoot, "interval");

    const uint &positionOutput = getSetting<uint>(cfgRoot, "positionOutput");
    const uint &positionInverval = getSetting<uint>(cfgRoot, "positionInterval");

    const uint &forceOutput = getSetting<uint>(cfgRoot, "forceOutput");
    const uint &dumpCoverage = getSetting<uint>(cfgRoot, "dumpCoverage");

    const uint nZerosBeforeTermination = 10000;

    rng.initialize(time(nullptr));

    const double xi = 1-exp(-1/ld);
    const double gamma0 = alpha*F0;

    SOSSolver solver(L, W, alpha, gamma0, true);

    setBoundariesFromIDs(&solver, {0,0,0,0}, L, W);

    RDLPotential rdlpotential(solver, s0, ld);
    solver.addLocalPotential(&rdlpotential);
    solver.registerObserver(&rdlpotential);

    ExtraNeighbor extraNeighbor(solver);
    solver.addLocalPotential(&extraNeighbor);

    const double Pl = ld*xi*F0;
    RDLExtraNeighborSurface rdlExtraSurface(solver, rdlpotential, extraNeighbor, Pl);

    ConfinedConstantConcentration diff(solver);

    Lattice lattice;

    lattice.addEvent(solver);
    lattice.addEvent(rdlExtraSurface);
    lattice.addEvent(diff);


    Coverage coverage(solver, dumpCoverage == 1, interval, nCycles);
    lattice.addEvent(coverage);

    Time time(solver);
    lattice.addEvent(time);

    DetectZeroCoverage detectZeroCoverage(&coverage, nZerosBeforeTermination);
    lattice.addEvent(detectZeroCoverage);

    DumpForceProfile forceProfile(rdlExtraSurface, interval);

    if (forceOutput == 1)
    {
        lattice.addEvent(forceProfile);
    }

    DumpSystem dumper(solver, positionInverval, path);

    if (positionOutput == 1)
    {
        lattice.addEvent(dumper);
    }

    H5Wrapper::Root h5root(path +  addProcEnding(argv, argc, "extraneighbor", "h5"));
    H5Wrapper::Member &simRoot = setuph5(h5root, getProc(argv, argc), L, W);

    StoreHeights storeHeights(solver, interval, simRoot);

    lattice.enableOutput(true, stdoutInterval);
    lattice.enableProgressReport();
    lattice.enableEventValueStorage(true,
                                    ignisOutput == 1,
                                    "ignisSOS.ign",
                                    "/tmp",
                                    interval);

    initializeSurface(solver, "random");
    const double h0 = rdlExtraSurface.getRdlEquilibrium();
    const int hm = solver.heights().max();

    if (h0 < hm + 1)
    {
        rdlExtraSurface.setHeight(hm+1);
    }
    else
    {
        rdlExtraSurface.setHeight(h0);
    }

    if (omegaSign == 0 && dumpCoverage == 1)
    {
        lattice.addEvent(storeHeights);
    }

    lattice.eventLoop(nCycles);


    /*
     *
     *
     *
     */

    simRoot["alpha"] = alpha;
    simRoot["omegaSign"] = omegaSign;
    simRoot["omegaVal"] = omegaVal;
    simRoot["F0"] = F0;
    simRoot["s0"] = s0;
    simRoot["ld"] = ld;
    simRoot["interval"] = interval;

    simRoot["eq_storedEventValues"] = lattice.storedEventValues();

    simRoot["eq_heights"] = solver.heights();

    if (dumpCoverage == 1)
    {
        simRoot["eq_coverage_matrix"] = coverage.coverage();
    }

    h5root.flush();

    if (omegaSign == 0)
    {
        return 0;
    }

    if (dumpCoverage == 1)
    {
        lattice.addEvent(storeHeights);
    }

    solver.setGamma(solver.gamma() + log(1 + omegaSign*omegaVal));

    diff.fixConcentration();

    lattice.eventLoop(nCyclesOmega);

    simRoot["omega_storedEventValues"] = lattice.storedEventValues();

    simRoot["omega_heighs"] = solver.heights();

    if (dumpCoverage)
    {
        simRoot["omega_coverage_matrix"] = coverage.coverage();
    }

    /*
     *
     *
     *
     */

    return 0;
}
