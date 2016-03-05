#include <SOSkMC.h>

#include "../apputils.h"

#include <HDF5Wrapper/include/hdf5wrapper.h>
#include <libconfig_utils/libconfig_utils.h>

class Cconvergence : public SOSEvent
{
public:

    Cconvergence(SOSSolver &solver,
                 FirstPassageContinuum &diffusionEvent,
                 const double delta,
                 const uint nThermCycles,
                 const uint nCyclesPerEquilibrium) :
        SOSEvent(solver, "cconv", "", true, true),
        m_mutexSolver(solver),
        m_diffusionEvent(diffusionEvent),
        m_a0(m_diffusionEvent.c()-delta),
        m_b0(m_diffusionEvent.c()+delta),
        m_nThermCycles(nThermCycles),
        m_nCyclesPerEquilibrium(nCyclesPerEquilibrium)
    {

    }

    void initialize()
    {
        m_count = 0;
        m_iteration = 1;

        m_a = m_a0;
        m_b = m_b0;

        m_faSet = false;
        m_fbSet = false;
    }

    void execute()
    {
        setValue(m_diffusionEvent.c());

        if (m_count == m_nThermCycles)
        {
            m_h0 = solver().averageHeight();
            m_T0 = solver().currentTime();
        }

        else if (m_count > m_nThermCycles)
        {
            m_v += (solver().averageHeight() - m_h0)/(solver().currentTime() - m_T0);
        }

        m_count++;
    }

    void reset()
    {
        BADAss(m_b, >, m_a);

        if (cycle() == 0)
        {
            return;
        }

        if (cycle() % m_nCyclesPerEquilibrium == 0)
        {
            m_v /= (m_count - m_nThermCycles - 1);

            if (m_b-m_a < 1E-3)
            {
                terminateLoop("root found.");
            }

            if (m_faSet)
            {
                if (m_v < m_fa)
                {
                    stringstream ss;
                    ss << "root found to noise precision: v="
                       << m_v << " va=" << m_fa << " c=" << m_diffusionEvent.c();

                    terminateLoop(ss.str());
                    return;
                }
            }

            if (m_fbSet)
            {
                if (m_v > m_fb)
                {
                    stringstream ss;
                    ss << "root found to noise precision: v="
                       << m_v << " vb=" << m_fb << " c=" << m_diffusionEvent.c();

                    terminateLoop(ss.str());
                    return;
                }
            }

            if (m_v < 0)
            {
                m_a = m_diffusionEvent.c();
                m_fa = m_v;
                m_faSet = true;
            }

            else
            {
                m_b = m_diffusionEvent.c();
                m_fb = m_v;
                m_fbSet = true;
            }

            m_diffusionEvent.setc((m_a+m_b)/2);

            m_count = 0;
            m_v = 0;

            m_iteration++;
        }
    }

    const double &a() const
    {
        return m_a;
    }

    const double &b() const
    {
        return m_b;
    }

    const double &fa() const
    {
        return m_fa;
    }

    const double &fb() const
    {
        return m_fb;
    }

    const bool &faSet() const
    {
        return m_faSet;
    }

    const bool &fbSet() const
    {
        return m_fbSet;
    }

    const double &v() const
    {
        return m_v;
    }

private:

    SOSSolver &m_mutexSolver;
    FirstPassageContinuum &m_diffusionEvent;

    double m_h0;
    double m_T0;
    uint m_count;

    const double m_a0;
    const double m_b0;

    double m_a;
    double m_b;
    double m_mid;

    double m_fa;
    double m_fb;
    double m_v;

    bool m_faSet;
    bool m_fbSet;

    uint m_iteration;
    const uint m_nThermCycles;
    const uint m_nCyclesPerEquilibrium;
};

int main(int argv, char **argc)
{
    rng.initialize(time(nullptr) % 1000000);
    //    rng.initialize(123456789);

    string cfgName = getCfgName(argv, argc, "cconv");

    Config cfg;
    cfg.readFile(cfgName.c_str());

    const Setting & cfgRoot = cfg.getRoot();

    const string path = getSetting<string>(cfgRoot, "path") + "/";
    const int type = getSetting<int>(cfgRoot, "type");

    const uint nSurfaceThermCycles = getSetting<uint>(cfgRoot, "nSurfaceThermCycles");
    const uint nThermCycles = getSetting<uint>(cfgRoot, "nThermCycles");
    const uint nCyclesPerEquilibrium = getSetting<uint>(cfgRoot, "nCyclesPerEquilibrium");
    const uint maxIterations = getSetting<uint>(cfgRoot, "maxIterations");
    const double c0 = getSetting<double>(cfgRoot, "c0");
    const double delta = getSetting<double>(cfgRoot, "delta");

    const int halfSize = getSetting<int>(cfgRoot, "halfSize");
    const double maxdt = getSetting<double>(cfgRoot, "maxdt");

    const uint L = getSetting<uint>(cfgRoot, "L");
    const uint W = getSetting<uint>(cfgRoot, "W");
    const double h0 = getSetting<double>(cfgRoot, "h0");
    const double alpha = getSetting<double>(cfgRoot, "alpha");


    const double gammaEq = 0;
    SOSSolver solver(L, W, alpha, gammaEq);

    setBoundariesFromIDs(&solver, {0,0,0,0}, L, W);

    ConstantConfinement confiningSurface(solver, h0);
    FirstPassageContinuum *diffusion;

    if (type == 0)
    {
        diffusion = new RadialFirstPassage(solver, maxdt, halfSize, c0);
    }
    else
    {
        diffusion = new AStarFirstPassage(solver, maxdt, halfSize, c0);
    }

    ParticleNumberConservator pnc(solver);

    Cconvergence converger(solver, *diffusion, delta, nThermCycles, nCyclesPerEquilibrium);

    Lattice lattice;
    //    lattice.enableOutput(false);
    lattice.enableProgressReport();
    lattice.enableEventValueStorage(false, true, "ignisSOS.ign", "/tmp", 100);

    lattice.addEvent(solver);
    lattice.addEvent(confiningSurface);
    lattice.addEvent(pnc);
    lattice.addEvent(diffusion);
    lattice.addEvent(converger);

    Time time(solver);
    ConcentrationTracker conc(solver);
    AverageHeight avgHeight(solver);

    lattice.addEvent(time);
    lattice.addEvent(conc);
    lattice.addEvent(avgHeight);

    initializeSurface(solver, "thermalized", 3 + type, nSurfaceThermCycles);

    lattice.eventLoop(nCyclesPerEquilibrium*maxIterations);

    H5Wrapper::Root h5root(path +  addProcEnding(argv, argc, "cconv", "h5"));

    stringstream sizeDesc;
    sizeDesc << L << "x" << W;
    H5Wrapper::Member &sizeRoot = h5root.addMember(sizeDesc.str());

    timeval tv;
    gettimeofday(&tv, nullptr);

    __int64_t run_ID = 1000*tv.tv_sec + tv.tv_usec/1000 + 10000000000000u*getProc(argv, argc);

    H5Wrapper::Member &simRoot = sizeRoot.addMember(run_ID);

    simRoot["type"] = type;
    simRoot["nCycles"] = solver.cycle();

    simRoot["L"] = L;
    simRoot["W"] = W;
    simRoot["alpha"] = alpha;
    simRoot["h0"] = h0;

    simRoot["nThermCycles"] = nThermCycles;
    simRoot["nCyclesPerEquilibrium"] = nCyclesPerEquilibrium;
    simRoot["c0"] = c0;
    simRoot["delta"] = delta;

    simRoot["halfSize"] = halfSize;
    simRoot["maxdt"] = maxdt;

    simRoot["c"] = diffusion->c();
    simRoot["a"] = converger.a();
    simRoot["b"] = converger.b();
    simRoot["fa"] = converger.fa();
    simRoot["fb"] = converger.fb();
    simRoot["faSet"] = converger.faSet();
    simRoot["fbSet"] = converger.fbSet();

    delete diffusion;

    return 0;
}
