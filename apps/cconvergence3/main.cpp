#include <SOSkMC.h>

#include "../apputils.h"

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
        m_delta(delta),
        m_nThermCycles(nThermCycles),
        m_nCyclesPerEquilibrium(nCyclesPerEquilibrium)
    {

    }

    void initialize()
    {
        m_count = 0;

        m_iteration = 1;
    }

    void execute()
    {
        setValue(m_diffusionEvent.c());

        if (m_count == m_nThermCycles)
        {
            m_h0 = solver().averageHeight();
            m_T0 = solver().currentTime();
        }

        m_count++;
    }

    void reset()
    {
        if (cycle() == 0)
        {
            return;
        }

        if (cycle() % m_nCyclesPerEquilibrium == 0)
        {
            const double v = (solver().averageHeight() - m_h0)/(solver().currentTime() - m_T0);
            int sign;

            if (m_iteration == 1)
            {
                m_v0 = v;
            }

            if (v > 0)
            {
                sign = -1;
            }

            else
            {
                sign = 1;
            }

            m_diffusionEvent.setc(m_diffusionEvent.c() + sign*m_delta/pow(m_iteration, 1.5));

            m_count = 0;

            m_iteration++;
        }
    }

private:

    SOSSolver &m_mutexSolver;
    FirstPassageContinuum &m_diffusionEvent;

    double m_h0;
    double m_T0;
    uint m_count;

    double m_v0;

    uint m_iteration;
    const double m_delta;

    const uint m_nThermCycles;
    const uint m_nCyclesPerEquilibrium;
};

int main(int argc, char **argv)
{
    double h0;
    double c0;

    if (argc == 3)
    {
        c0 = atof(argv[1]);
        h0 = atof(argv[2]);
    }

    else
    {
        c0 = 0.1;
        h0 = 10.;
    }

    rng.initialize(time(nullptr) % 1000000);

    const uint L = 30;
    const uint W = 30;
    const double alpha = 1.0;
    const double gammaEq = 0;
    SOSSolver solver(L, W, alpha, gammaEq);

    setBoundariesFromIDs(&solver, {0,0,0,0}, L, W);

    ConstantConfinement confiningSurface(solver, h0);

    const double maxdt = 0.01;
    RadialFirstPassage diffusion(solver, maxdt, 3, c0);

    ParticleNumberConservator pnc(solver);

    const uint nThermCycles = 2000;
    const uint nCyclesPerEquilibrium = 5000;
    const double delta = 0.05;
    Cconvergence converger(solver, diffusion, delta, nThermCycles, nCyclesPerEquilibrium);

    Lattice lattice;
    //    lattice.enableOutput(false);
    lattice.enableProgressReport();
    lattice.enableEventValueStorage(false, true, "ignisSOS.ign", "/tmp", 1);

    lattice.addEvent(solver);
    lattice.addEvent(confiningSurface);
    lattice.addEvent(diffusion);
    lattice.addEvent(converger);
    lattice.addEvent(pnc);

    Time time(solver);
    ConcentrationTracker conc(solver);
    AverageHeight avgHeight(solver);

    lattice.addEvent(time);
    lattice.addEvent(conc);
    lattice.addEvent(avgHeight);

    initializeSurface(solver, "thermalized", 3, 3000);

    const uint nIterations = 20;
    lattice.eventLoop(nCyclesPerEquilibrium*nIterations);

    return 0;
}
