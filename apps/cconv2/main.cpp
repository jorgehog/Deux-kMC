#include <SOSkMC.h>

#include "../apputils.h"

#include <HDF5Wrapper/include/hdf5wrapper.h>
#include <libconfig_utils/libconfig_utils.h>


class MeasureSomething : public SOSEvent
{
public:
    MeasureSomething(const SOSSolver &solver,
                     const std::function<double(const SOSSolver &)> f) :
        SOSEvent(solver, "Something", "", true),
        m_f(f)
    {

    }

    void execute()
    {
        this->setValue(m_f(solver()));
    }

private:
    const std::function<double(const SOSSolver &)> &m_f;
};

class AvgLast : public SOSEvent
{
public:

    AvgLast(const SOSEvent &event,
            const uint N,
            bool hasOutput,
            bool storeValue) :
        SOSEvent(event.solver(),
                 event.type() + "Avg",
                 event.unit(),
                 hasOutput,
                 storeValue),
        m_event(event),
        m_N(N),
        m_values(N)
    {
        setDependency(event);
    }

private:

    const SOSEvent &m_event;

    const uint m_N;
    vec m_values;

    uint m_count;
    double m_sum;

    // Event interface
public:
    void execute()
    {
        const double value = m_event.value();

        if (m_count < m_N)
        {
            m_values(m_count) = value;
            m_sum += value;
        }

        else
        {
            for (uint i = 0; i < m_N - 1; ++i)
            {
                m_values(i+1) = m_values(i);
            }

            m_sum += value - m_values(0);
            m_values(0) = value;

            setValue(m_sum/m_N);
        }

        m_count++;

    }

    void initialize()
    {
        m_count = 0;
        m_sum = 0;
    }

};


//source:
//http://www.johndcook.com/blog/2012/06/14/root-finding-with-noisy-functions/
class BisectC
{
public:
    BisectC(const uint L,
            const uint W,
            const double alpha,
            const double h0,
            const int type,
            const uint nThermCycles,
            const uint nCyclesPerEquilibrium,
            const bool toFile) :
        m_L(L),
        m_W(W),
        m_alpha(alpha),
        m_h0(h0),
        m_type(type),
        m_nThermCycles(nThermCycles),
        m_nCyclesPerEquilibrium(nCyclesPerEquilibrium),
        m_toFile(toFile)
    {

    }

    void bisect(double &a,
                double &b,
                double fa,
                double fb,
                const double tol = 1e-3)
    {
        if (b - a < tol)
        {
            cout << "tolerance reached: " << a << " " << b << endl;
            return;
        }

        const double mid = (b+a)/2;
        const double fmid = estimateSpeed(mid);

        if (fmid < fa || fmid > fb)
        {
            cout << "error range reached: " << fa << " " << fmid << " " << fb << endl;
            return;
        }

        if (fmid < 0)
        {
            a = mid;
            fa = fmid;
        }

        else
        {
            b = mid;
            fb = fmid;
        }

        return bisect(a, b, fa, fb, tol);
    }

    void bisect(double &a,
                double &b,
                const double tol=1e-3)
    {
        const double fa = estimateSpeed(a);
        const double fb = estimateSpeed(b);

        if (fa*fb >= 0)
        {
            throw std::runtime_error("zero is not in interval.");
        }

        bisect(a, b, fa, fb, tol);
    }

    double estimateSpeed(const double c)
    {
        cout << "running " << c << endl;

        SOSSolver solver(m_L, m_W, m_alpha, 0);

        vector<vector<Boundary *> > b = setBoundariesFromIDs(&solver, {0,0,0,0}, m_L, m_W);

        ConstantConfinement confiningSurface(solver, m_h0);
        FirstPassageContinuum *diffusion;

        const double maxdt = 0.01;
        const int halfSize = 3;

        if (m_type == 0)
        {
            diffusion = new RadialFirstPassage(solver, maxdt, halfSize, c);
        }
        else
        {
            diffusion = new AStarFirstPassage(solver, maxdt, halfSize, c);
        }


        ParticleNumberConservator pnc(solver);

        diffusion->setDependency(pnc);

        GrowthSpeed speed(solver);
        speed.setOnsetTime(m_nThermCycles);
        AverageHeight avgHeight(solver);


        Lattice lattice;

        lattice.enableOutput(true, m_nCyclesPerEquilibrium/10);
//        lattice.enableOutput(true, 1);
        lattice.enableProgressReport(true);
        lattice.enableEventValueStorage(false, m_toFile, "ignisSOS.ign", "/tmp", 100);

        lattice.addEvent(solver);
        lattice.addEvent(confiningSurface);
        lattice.addEvent(pnc);
        lattice.addEvent(diffusion);
        lattice.addEvent(avgHeight);
        lattice.addEvent(speed);

//        lattice.addEvent(new MeasureSomething(solver, [&] (const SOSSolver &s)
//        {
//            return s.confiningSurfaceEvent().height() - s.averageHeight();
//        }));

//        DumpSystem dump(solver, 1);
//        lattice.addEvent(dump);

        initializeSurface(solver, "random", 3 + m_type, m_nThermCycles, false);

        lattice.eventLoop(m_nCyclesPerEquilibrium);

        const double v = speed.value();

        delete diffusion;

        delete b[0][0];
        delete b[0][1];
        delete b[1][0];
        delete b[1][1];

        return v;
    }

private:

    const uint m_L;
    const uint m_W;
    const double m_alpha;

    const double m_h0;
    const int m_type;

    const uint m_nThermCycles;
    const uint m_nCyclesPerEquilibrium;

    const bool m_toFile;

};

int main(int argv, char **argc)
{
    rng.initialize(time(nullptr) % 1000000);
//            rng.initialize(361632621);

    string cfgName = getCfgName(argv, argc, "cconv2");

    Config cfg;
    cfg.readFile(cfgName.c_str());

    const Setting & cfgRoot = cfg.getRoot();

    const string path = getSetting<string>(cfgRoot, "path") + "/";
    const int type = getSetting<int>(cfgRoot, "type");

    const uint nThermCycles = getSetting<uint>(cfgRoot, "nThermCycles");
    const uint nCyclesPerEquilibrium = getSetting<uint>(cfgRoot, "nCyclesPerEquilibrium");

    const double a0 = getSetting<double>(cfgRoot, "a");
    const double b0 = getSetting<double>(cfgRoot, "b");

    //    const int halfSize = getSetting<int>(cfgRoot, "halfSize");
    //    const double maxdt = getSetting<double>(cfgRoot, "maxdt");

    const uint L = getSetting<uint>(cfgRoot, "L");
    const uint W = getSetting<uint>(cfgRoot, "W");
    const double h0 = getSetting<double>(cfgRoot, "h0");
    const double alpha = getSetting<double>(cfgRoot, "alpha");

    const bool toFile = getSetting<uint>(cfgRoot, "toFile") == 1;

    BisectC bc(L, W, alpha, h0, type,
               nThermCycles, nCyclesPerEquilibrium,
               toFile);

    double a = a0;
    double b = b0;

    bc.bisect(a, b);

    cout << (a + b)/2 << endl;

    H5Wrapper::Root h5root(path +  addProcEnding(argv, argc, "cconv2", "h5"));

    stringstream sizeDesc;
    sizeDesc << L << "x" << W;
    H5Wrapper::Member &sizeRoot = h5root.addMember(sizeDesc.str());

    timeval tv;
    gettimeofday(&tv, nullptr);

    __int64_t run_ID = 1000*tv.tv_sec + tv.tv_usec/1000 + 10000000000000u*getProc(argv, argc);

    H5Wrapper::Member &simRoot = sizeRoot.addMember(run_ID);

    simRoot["type"] = type;

    simRoot["L"] = L;
    simRoot["W"] = W;
    simRoot["alpha"] = alpha;
    simRoot["h0"] = h0;

    simRoot["nThermCycles"] = nThermCycles;
    simRoot["nCyclesPerEquilibrium"] = nCyclesPerEquilibrium;

    //    simRoot["halfSize"] = halfSize;
    //    simRoot["maxdt"] = maxdt;

    simRoot["a"] = a;
    simRoot["b"] = b;

    return 0;
}
