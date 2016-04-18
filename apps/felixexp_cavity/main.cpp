#include <SOSkMC.h>

#include "../apputils.h"

using ignis::Lattice;

class PartialNeighbors : public LocalPotential, public Observer<Subjects>
{
public:
    PartialNeighbors(SOSSolver &solver,
                     AverageHeightBoundary &ahb) :
        LocalPotential(solver),
        m_ahb(ahb)
    {

    }

    virtual ~PartialNeighbors() {}

private:

    AverageHeightBoundary &m_ahb;

    double m_prevAverage;

    // LocalPotential interface
public:

    static constexpr double Eb = 1.0;

    double potential(const uint x, const uint y) const
    {
        (void) y;

        const bool atBoundary = x == 0;

        if (atBoundary)
        {
            const int hBoundary = round(m_ahb.average());
            const double &hAvg = m_ahb.average();

            if (hAvg < hBoundary)
            {
                return 0;
            }

            else
            {
                const double probFilled = (hAvg - hBoundary);
                const double corr = probFilled*Eb;

                //This means that there is a contribution Eb already
                if (probFilled >= 0.5)
                {
                    BADAssBool(m_ahb.isBlockedLattice(-1, y, solver().height(x, y)), "hmm", [&] ()
                    {
                        BADAssSimpleDump(corr, m_ahb.average());
                    });

                    return corr - Eb;
                }

                else
                {
                    return corr;
                }
            }
        }

        return 0;
    }

    // Observer interface
public:
    void initializeObserver(const Subjects &subject)
    {
        (void) subject;

        m_prevAverage = m_ahb.average();
    }
    void notifyObserver(const Subjects &subject)
    {
        (void) subject;

        if (m_ahb.average() != m_prevAverage)
        {
            m_ahb.affectSurfaceSites();
            m_prevAverage = m_ahb.average();
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

        if (xi >= solver().length() - 0.5)
        {

            return 2*solver().length() - xi - 1;
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


class ReflAvgHybrid : public AverageHeightBoundary
{
public:

    ReflAvgHybrid(SOSSolver &solver, const uint averageHeightDepth) :
        AverageHeightBoundary(solver, averageHeightDepth, 0,
                              solver.length(), solver.width(),
                              Boundary::orientations::LAST, solver.length() - 1)
    {

    }

    // Boundary interface
public:
    double transformContinousCoordinate(const double xi, const double xj, const double xk) const
    {
        (void) xj;
        (void) xk;

        if (xi >= solver().length() - 0.5)
        {

            return 2*solver().length() - xi - 1;
        }

        else
        {
            return xi;
        }
    }

    bool isBlockedContinous(const double xi, const double xj, const double xk) const
    {
        (void) xi;
        (void) xj;
        (void) xk;

        return false;
    }

    void closestImage(const double xi, const double xj, const double xk, const double xti, const double xtj, const double xtk, double &dxi, double &dxj, double &dxk) const
    {
        noImage(xi, xj, xk, xti, xtj, xtk, dxi, dxj, dxk);
    }
};


class StoreHeights : public SOSEvent
{
public:
    StoreHeights(const SOSSolver &solver, const uint interval, H5Wrapper::Member &h5group) :
        SOSEvent(solver, "storeheights"),
        m_interval(interval),
        m_h5group(h5group.addMember("stored_heights"))
    {

    }

private:
    const uint m_interval;
    H5Wrapper::Member &m_h5group;

    // Event interface
public:
    void execute()
    {
        if (cycle() % m_interval == 0)
        {
            m_h5group[to_string(cycle())] = solver().heights();
        }
    }
};

class StoreParticles : public SOSEvent
{
public:
    StoreParticles(const SOSSolver &solver,
                 const OfflatticeMonteCarlo &omc,
                 const uint interval,
                 H5Wrapper::Member &h5group) :
        SOSEvent(solver, "storeparticles"),
        m_omc(omc),
        m_interval(interval),
        m_h5group(h5group.addMember("stored_particles"))
    {

    }

private:
    const OfflatticeMonteCarlo &m_omc;
    const uint m_interval;
    H5Wrapper::Member &m_h5group;

    // Event interface
public:
    void execute()
    {
        if (cycle() % m_interval == 0)
        {
            m_h5group[to_string(cycle())] = m_omc.particlePositions()(span::all,
                                                                      span(0, m_omc.nOfflatticeParticles() - 1)).eval();
        }
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

    const uint &nCycles = getSetting<uint>(cfgRoot, "nCycles");
    const uint &interval = getSetting<uint>(cfgRoot, "interval");
    const uint &output = getSetting<uint>(cfgRoot, "output");

    rng.initialize(time(nullptr));
//        rng.initialize(100101010);

    const string initType = "thermalized";
//        const string initType = "random";

    const double gamma = log(1 + omega);

    SOSSolver solver(L, W, alpha, gamma, true);

    AverageHeightBoundary x0(solver, 5, 0, L, W, Boundary::orientations::FIRST, 0);
    ReflAvgHybrid x1(solver, 5);

    Periodic y0(W, Boundary::orientations::FIRST);
    Periodic y1(W, Boundary::orientations::LAST);

    solver.setBoundaries({{&x0, &x1}, {&y0, &y1}});

    RadialFirstPassage diff(solver, 0.01, 3, getMFPTConstant(height, alpha, 0));
    FixedSurface confSurface(solver, height);

    solver.addConcentrationBoundary(0, Boundary::orientations::FIRST, omega);

    PartialNeighbors pn(solver, x0);
    solver.addLocalPotential(&pn);
    solver.registerObserver(&pn);

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

    if (output)
    {
        lattice.addEvent(dumper);
    }
    else
    {
        lattice.addEvent(storeHeights);
        lattice.addEvent(storeParticles);
    }

    initializeSurface(solver, initType, 3, 10000);

    lattice.eventLoop(nCycles);

    simRoot["alpha"] = alpha;
    simRoot["omega"] = omega;
    simRoot["height"] = height;
    simRoot["time"] = colvec(lattice.storedEventValues().col(0));

    return 0;
}
