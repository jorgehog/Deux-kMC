#include <SOSkMC.h>

#include "../apputils.h"

#include "extraneighbor.h"
#include "rdlextraneighborsurface.h"

using ignis::Lattice;

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

    const double Pl = getSetting<double>(cfgRoot, "Pl");
    const double alpha = getSetting<double>(cfgRoot, "alpha");
    const double omegaShift = getSetting<double>(cfgRoot, "omegaShift");

    const double ld = getSetting<double>(cfgRoot, "ld");
    const double s0 = getSetting<double>(cfgRoot, "s0");

    const uint nCycles = getSetting<uint>(cfgRoot, "nCycles");
    const uint interval = getSetting<uint>(cfgRoot, "interval");
    const uint ignisOutput = getSetting<uint>(cfgRoot, "ignisOutput");
    const uint positionOutput = getSetting<uint>(cfgRoot, "positionOutput");
    const uint dumpCoverage = getSetting<uint>(cfgRoot, "dumpCoverage");

    const uint nZerosBeforeTermination = 10000;

    rng.initialize(time(nullptr));

    const double eTerm = 1-exp(-1/ld);
    const double gamma0 = alpha*Pl/(ld*eTerm);

    SOSSolver solver(L, W, alpha, gamma0, true);

    setBoundariesFromIDs(&solver, {0,0,0,0}, L, W);

    RDLPotential rdlpotential(solver, s0, ld);
    solver.addLocalPotential(&rdlpotential);
    solver.registerObserver(&rdlpotential);

    ExtraNeighbor extraNeighbor(solver);
    solver.addLocalPotential(&extraNeighbor);

    RDLExtraNeighborSurface rdlSurface(solver, rdlpotential, extraNeighbor, Pl);

    ConfinedConstantConcentration diff(solver);


    Lattice lattice;

    lattice.addEvent(solver);
    lattice.addEvent(rdlSurface);
    lattice.addEvent(diff);

    Coverage coverage(solver, dumpCoverage == 1, interval, nCycles);
    lattice.addEvent(coverage);

    Time time(solver);
    lattice.addEvent(time);

    DetectZeroCoverage detectZeroCoverage(&coverage, nZerosBeforeTermination);
    lattice.addEvent(detectZeroCoverage);

    DumpSystem dumper(solver, interval, path);

    if (positionOutput)
    {
        lattice.addEvent(dumper);
    }

    lattice.enableOutput(true, interval);
    lattice.enableProgressReport();
    lattice.enableEventValueStorage(true,
                                    ignisOutput,
                                    "ignisSOS.ign",
                                    "/tmp",
                                    interval);

    initializeSurface(solver, "random");
    const double h0 = rdlSurface.getRdlEquilibrium();
    const int hm = solver.heights().max();

    if (h0 < hm + 1)
    {
        rdlSurface.setHeight(hm+1);
    }
    else
    {
        rdlSurface.setHeight(h0);
    }

    lattice.eventLoop(nCycles);


    /*
     *
     *
     *
     */

    H5Wrapper::Root h5root(path +  addProcEnding(argv, argc, "extraneighbor", "h5"));
    H5Wrapper::Member &simRoot = setuph5(h5root, getProc(argv, argc), L, W);

    simRoot["alpha"] = alpha;
    simRoot["omegaShift"] = omegaShift;
    simRoot["Pl"] = Pl;
    simRoot["s0"] = s0;
    simRoot["interval"] = interval;

    simRoot["eq_storedEventValues"] = lattice.storedEventValues();

    if (dumpCoverage == 1)
    {
        simRoot["eq_coverage_matrix"] = coverage.coverage();
    }

    if (omegaShift == 0)
    {
        return 0;
    }

    solver.setGamma(solver.gamma() + log(1 + omegaShift));

    lattice.eventLoop(nCycles);

    simRoot["omega_storedEventValues"] = lattice.storedEventValues();

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
