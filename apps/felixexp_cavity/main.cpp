#include <SOSkMC.h>

#include "../apputils.h"

using ignis::Lattice;

class InFluxer : public SOSEvent
{
public:
    InFluxer(const SOSSolver &solver,
             OfflatticeMonteCarlo &omc,
             const uint influxInterval) :
        SOSEvent(solver, "InFluxer", "", true, true),
        m_x0(0),
        m_omc(omc),
        m_influxInterval(influxInterval),
        m_boundaryHeights(solver.width()),
        m_accuBoundaryHeights(solver.width())
    {

    }

private:
    const uint m_x0;
    OfflatticeMonteCarlo &m_omc;

    const uint m_influxInterval;

    vector<double> m_boundaryHeights;
    vector<double> m_accuBoundaryHeights;

    // Event interface
public:
    void execute()
    {
        if (cycle() == 0)
        {
            return;
        }

        if (cycle() % m_influxInterval == 0)
        {
            double dhSum = 0;
            const double &hl = solver().confiningSurfaceEvent().height();
            for (uint y = 0; y < solver().width(); ++y)
            {
                const int &h = solver().height(m_x0, y);

                double dh = hl - h - 2;
                m_boundaryHeights.at(y) = dh;

                dhSum += dh;
                m_accuBoundaryHeights.at(y) = dhSum;
            }

            if (dhSum <= 0)
            {
                terminateLoop("Inlet clogged. Exiting");
                return;
            }

            const double R = rng.uniform()*dhSum;

            const uint y0 = binarySearchAndScan(&m_accuBoundaryHeights.front(),
                                                solver().width(),
                                                R);

            const double z0 = solver().height(m_x0, y0) + 1 + m_accuBoundaryHeights.at(y0) - R;

            m_omc.insertParticle(m_x0, y0, z0);
        }
    }
};

class ReflConstantHybrid : public SOSBoundary
{
public:

    ReflConstantHybrid(SOSSolver &solver) :
        SOSBoundary(solver, Boundary::orientations::LAST)
    {

    }

    // Boundary interface
public:
    double transformContinousCoordinate(const double xi, const double xj, const double xk) const
    {
        (void) xj;
        (void) xk;

        if (xi >= solver().length())
        {
            return 2*(solver().length() - 1) - xi + 1;
        }

        else
        {
            return xi;
        }
    }

    int transformLatticeCoordinate(const int xi, const int xj, const int xk) const
    {
        (void) xj;
        (void) xk;

        return xi;
    }

    bool isBlockedContinous(const double xi, const double xj, const double xk) const
    {
        (void) xi;
        (void) xj;
        (void) xk;

        return false;
    }

    bool isBlockedLattice(const int xi, const int xj, const int xk) const
    {
        (void) xj;

        const int loc = solver().length();
        const int h = 0;

        return (xi >= loc) && (xk <= h);
    }

    void closestImage(const double xi, const double xj, const double xk, const double xti, const double xtj, const double xtk, double &dxi, double &dxj, double &dxk) const
    {
        noImage(xi, xj, xk, xti, xtj, xtk, dxi, dxj, dxk);
    }
};


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

    const uint &influxInterval = getSetting<uint>(cfgRoot, "influxInterval");

    const uint &nCycles = getSetting<uint>(cfgRoot, "nCycles");
    const uint &interval = getSetting<uint>(cfgRoot, "interval");
    const uint &output = getSetting<uint>(cfgRoot, "output");

    rng.initialize(time(nullptr));
//    rng.initialize(100101010);

    const double gamma0 = 0;
    const double gamma = log(1 + omega) + gamma0;

    SOSSolver solver(L, W, alpha, gamma, true);

    Boundary *x0;
    Boundary *x1;

    x0 = new AverageHeightBoundary(solver, 5, 0, L, W, Boundary::orientations::FIRST, 0);
    x1 = new ReflConstantHybrid(solver);

    Periodic y0(W, Boundary::orientations::FIRST);
    Periodic y1(W, Boundary::orientations::LAST);

    solver.setBoundaries({{x0, x1}, {&y0, &y1}});

    RadialFirstPassage diff(solver, 0.01, 3, getMFPTConstant(height, alpha, 0));
    FixedSurface confSurface(solver, height);

    Lattice lattice;

    lattice.addEvent(solver);
    lattice.addEvent(confSurface);
    lattice.addEvent(diff);

    InFluxer inFluxer(solver, diff, influxInterval);
    lattice.addEvent(inFluxer);

    DumpSystem dumper(solver, 100, path);
    lattice.addEvent(dumper);

    AverageHeight avgH(solver);
    lattice.addEvent(avgH);

    lattice.enableOutput(true, interval);
    lattice.enableProgressReport();
    lattice.enableEventValueStorage(true,
                                    output,
                                    "ignisSOS.ign",
                                    "/tmp",
                                    interval);

    //    const string initType = "thermalized";
    const string initType = "random";
    initializeSurface(solver, initType, 3, 10000);

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
    simRoot["omega"] = omega;
    simRoot["height"] = height;
    simRoot["coverage"] = colvec(lattice.storedEventValues().col(0));

    return 0;
}
