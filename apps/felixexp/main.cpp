#include <SOSkMC.h>

#include "../apputils.h"

#include <HDF5Wrapper/include/hdf5wrapper.h>
#include <libconfig_utils/libconfig_utils.h>


class InsertStep : public LatticeEvent
{
public:
    InsertStep(SOSSolver &solver, const uint interval, const uint size) :
        LatticeEvent("InsertStep"),
        m_solver(solver),
        m_interval(interval),
        m_size(size)
    {

    }

private:

    SOSSolver &m_solver;
    const uint m_interval;
    const uint m_size;
    uint m_currentLevel;

    void insertStep()
    {
        m_currentLevel += 1;

        for (uint x = solver().length()-m_size; x < solver().length(); ++x)
        {
            for (uint y = 0; y < m_size; ++y)
            {
                m_solver.setHeight(x, y, m_currentLevel, false);
            }
        }

        m_solver.calculateHeightDependentValues();
    }

    SOSSolver &solver() const
    {
        return m_solver;
    }

    // Event interface
public:
    void execute()
    {

    }

    void initialize()
    {
        m_currentLevel = 0;
        insertStep();
    }

    void reset()
    {
        if (cycle() == 0)
        {
            return;
        }

        if (cycle() % m_interval == 0)
        {
            insertStep();
        }
    }
};


int main(int argv, char **argc)
{
    rng.initialize(time(nullptr));
    //    rng.initialize(139444);

    string cfgName = getCfgName(argv, argc, "felixexp");

    Config cfg;
    cfg.readFile(cfgName.c_str());

    const Setting & cfgRoot = cfg.getRoot();

    const string path = getSetting<string>(cfgRoot, "path") + "/";

    const int halfSize = getSetting<int>(cfgRoot, "halfSize");
    const double maxdt = getSetting<double>(cfgRoot, "maxdt");

    const uint L = getSetting<uint>(cfgRoot, "L");
    const uint W = getSetting<uint>(cfgRoot, "W");
    const double h0 = getSetting<double>(cfgRoot, "h0");
    const double alpha = getSetting<double>(cfgRoot, "alpha");

    const uint size = getSetting<uint>(cfgRoot, "size");
    const uint interval = getSetting<uint>(cfgRoot, "interval");

    const uint boundaryDepth = getSetting<uint>(cfgRoot, "boundaryDepth");

    const uint nCycles = getSetting<uint>(cfgRoot, "nCycles");


    SOSSolver solver(L, W, alpha, 0);

    setBoundariesFromIDs(&solver, {5, 5, 5, 5}, L, W, 0, boundaryDepth);

    const double cBath = 2.0;
    ConcentrationProfile concProfile(solver, [&L, &cBath] (const uint x, const uint y)
    {
        return cBath - (cBath - 1)*x/(L-1);
    });

    ConstantConfinement confinement(solver, h0);

    InsertStep insertStep(solver, interval, size);

    Lattice lattice;
    lattice.addEvent(solver);
    lattice.addEvent(insertStep);
    lattice.addEvent(confinement);
    lattice.addEvent(concProfile);

    DumpSystem dumper(solver, 1, path);

    lattice.addEvent(dumper);


    //
    lattice.eventLoop(nCycles);
    //


    H5Wrapper::Root h5root(path +  addProcEnding(argv, argc, "felixexp", "h5"));

    stringstream sizeDesc;
    sizeDesc << L << "x" << W;
    H5Wrapper::Member &sizeRoot = h5root.addMember(sizeDesc.str());

    timeval tv;
    gettimeofday(&tv, nullptr);

    __int64_t run_ID = 1000*tv.tv_sec + tv.tv_usec/1000 + 10000000000000u*getProc(argv, argc);

    H5Wrapper::Member &simRoot = sizeRoot.addMember(run_ID);

    simRoot["L"] = L;
    simRoot["W"] = W;
    simRoot["alpha"] = alpha;
    simRoot["h0"] = h0;

    simRoot["halfSize"] = halfSize;
    simRoot["maxdt"] = maxdt;

    return 0;
}
