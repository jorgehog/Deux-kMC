#include <SOSkMC.h>

#include "../apputils.h"

#include "../extraneighbor/extraneighbor.h"

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
        bool atBoundary;

        if (m_ahb.dim() == 0)
        {
            atBoundary = x == m_ahb.location();
        }

        else
        {
            atBoundary = y == m_ahb.location();
        }

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

class VS : public SOSEvent
{
public:
    VS(const SOSSolver &solver) :
        SOSEvent(solver, "VS", "", true)
    {

    }

    // Event interface
public:
    void execute()
    {
        setValue(solver().volume() - solver().diffusionEvent().numberOfParticles());
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
    const double &omegaEq = getSetting<double>(cfgRoot, "omegaEq");
    const double &height = getSetting<double>(cfgRoot, "height");

    const uint &boundaryDepth = getSetting<uint>(cfgRoot, "boundaryDepth");

    const uint &nCycles = getSetting<uint>(cfgRoot, "nCycles");
    const uint &interval = getSetting<uint>(cfgRoot, "interval");
    const uint &output = getSetting<uint>(cfgRoot, "output");

    rng.initialize(time(nullptr));
    const string initType = "flat";
//    const string initType = "random";

    //        rng.initialize(100101010);

    const double gamma = log(1 + omega);

    SOSSolver solver(L, W, alpha, gamma, true);

    AverageHeightBoundary x0(solver, boundaryDepth, 0, L, W, Boundary::orientations::FIRST, 0);
    AverageHeightBoundary x1(solver, boundaryDepth, 0, L, W, Boundary::orientations::LAST, L-1);

    Periodic y0(W, Boundary::orientations::FIRST);
    Periodic y1(W, Boundary::orientations::LAST);

    solver.setBoundaries({{&x0, &x1}, {&y0, &y1}});

    const uint depositionBoxHalfSize = 3;
    const double maxDt = 0.01;

    RadialFirstPassage diff(solver, maxDt, depositionBoxHalfSize, getMFPTConstant(height, alpha, 0));
    FixedSurface confSurface(solver, height);

    solver.addConcentrationBoundary(0, Boundary::orientations::FIRST, omega);
    solver.addConcentrationBoundary(0, Boundary::orientations::LAST, omegaEq);

    PartialNeighbors pnLeft(solver, x0);
    solver.addLocalPotential(&pnLeft);
    solver.registerObserver(&pnLeft);

    PartialNeighbors pnRight(solver, x1);
    solver.addLocalPotential(&pnRight);
    solver.registerObserver(&pnRight);

    ExtraNeighbor en(solver);
    solver.addLocalPotential(&en);

    Time time(solver);
    VS vs(solver);

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
    lattice.addEvent(vs);

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
