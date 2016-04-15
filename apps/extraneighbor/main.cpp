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

int main(int argv, char** argc)
{
    string cfgName = getCfgName(argv, argc, "extraneighbor");

    Config cfg;
    cfg.readFile(cfgName.c_str());

    const Setting & cfgRoot = cfg.getRoot();

    const string path = getSetting<string>(cfgRoot, "path") + "/";
    const double Pl = getSetting<double>(cfgRoot, "Pl");
    const double alpha = getSetting<double>(cfgRoot, "alpha");
    const double omega = getSetting<double>(cfgRoot, "omega");

    const double ld = getSetting<double>(cfgRoot, "ld");
    const double s0 = getSetting<double>(cfgRoot, "s0");

    const uint nCycles = getSetting<uint>(cfgRoot, "nCycles");
    const uint interval = getSetting<uint>(cfgRoot, "interval");
    const uint output = getSetting<uint>(cfgRoot, "output");
    const uint dumpCoverage = getSetting<uint>(cfgRoot, "dumpCoverage");

    rng.initialize(time(nullptr));

    const uint L = 30;
    const uint W = 30;

    const double eTerm = 1-exp(-1/ld);
    const double gamma0 = alpha*Pl/(ld*eTerm);

    SOSSolver solver(L, W, alpha, gamma0, true);

    setBoundariesFromIDs(&solver, {0,0,0,0}, L, W);

    RDLPotential rdlpotential(solver, s0, ld);
    solver.addLocalPotential(&rdlpotential);
    solver.registerObserver(&rdlpotential);

    ExtraNeighbor extraNeighbor(solver, s0/eTerm);
    solver.addLocalPotential(&extraNeighbor);

    RDLExtraNeighborSurface rdlSurface(solver, rdlpotential, extraNeighbor, Pl);

    Diffusion *diff;

    diff = new ConfinedConstantConcentration(solver);

    Lattice lattice;

    lattice.addEvent(solver);
    lattice.addEvent(rdlSurface);
    lattice.addEvent(diff);

    Coverage coverage(solver, (dumpCoverage == 1) && (omega == 0), interval, nCycles);
    lattice.addEvent(coverage);

    DumpSystem dumper(solver, interval, path);

    if (output)
    {
        lattice.addEvent(dumper);
    }

    lattice.enableOutput(true, interval);
    lattice.enableProgressReport();
    lattice.enableEventValueStorage(true,
                                    output,
                                    "ignisSOS.ign",
                                    "/tmp",
                                    interval);

    initializeSurface(solver, "random");
    const double h0 = rdlSurface.getRdlEquilibrium();
    const int hm = solver.heights().max();

    if (hm > h0 - 1)
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

    double gamma2 = solver.gamma() + log(1 + omega);
    SOSSolver solver2(L, W, alpha, gamma2, true);

    setBoundariesFromIDs(&solver2, {0,0,0,0}, L, W);

    RDLPotential rdlpotential2(solver2, s0, ld);
    solver2.addLocalPotential(&rdlpotential2);
    solver2.registerObserver(&rdlpotential2);

    ExtraNeighbor extraNeighbor2(solver2, s0/eTerm);
    solver2.addLocalPotential(&extraNeighbor2);

    RDLExtraNeighborSurface rdlSurface2(solver2, rdlpotential2, extraNeighbor2, Pl);

    Diffusion *diff2;

    diff2 = new ConstantConcentration(solver2);

    Lattice lattice2;

    lattice2.addEvent(solver2);
    lattice2.addEvent(rdlSurface2);
    lattice2.addEvent(diff2);

    Coverage coverage2(solver2, dumpCoverage == 1, interval, nCycles);
    lattice2.addEvent(coverage2);

    DumpSystem dumper2(solver2, interval, path);

    if (output)
    {
        lattice2.addEvent(dumper2);
    }

    lattice2.enableOutput(true, interval);
    lattice2.enableProgressReport();
    lattice2.enableEventValueStorage(true,
                                     output,
                                     "ignisSOS2.ign",
                                     "/tmp",
                                     interval);

    solver2.setHeights(solver.heights(), false);
    rdlSurface2.setHeight(rdlSurface.height());

    double h;
    vec cov;

    if (omega == 0)
    {
        h = rdlSurface.height();
        cov = lattice.storedEventValues().col(0);
    }

    else
    {
        lattice2.eventLoop(nCycles);
        h = rdlSurface2.height();
        cov = lattice2.storedEventValues().col(0);
    }


    /*
     *
     *
     *
     */


    H5Wrapper::Root h5root(path +  addProcEnding(argv, argc, "extraneighbor", "h5"));
    H5Wrapper::Member &simRoot = setuph5(h5root, getProc(argv, argc), L, W);

    simRoot["alpha"] = alpha;
    simRoot["omega"] = omega;
    simRoot["Pl"] = Pl;
    simRoot["s0"] = s0;
    simRoot["h"] = h;
    simRoot["coverage"] = cov;
    simRoot["interval"] = interval;

    if (dumpCoverage == 1)
    {
        simRoot["coverage_matrix_eq"] = coverage.coverage();

        if (omega != 0)
        {
            simRoot["coverage_matrix"] = coverage2.coverage();
        }
    }

    return 0;
}
