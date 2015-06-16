#include "eqmu.h"

#include "solidonsolidreaction.h"

#include "miscevents.h"


void EqMu::initialize()
{
    resetCounters();

    update();
}

void EqMu::execute()
{
    setValue(dMu());
}

void EqMu::reset()
{
    update();
}

void EqMu::restart()
{
    resetCounters();
}

void EqMu::update()
{
    double localDissolutionRate = 0;

    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            localDissolutionRate += solver().reaction(x, y).dissolutionRate();
        }
    }

    localDissolutionRate /= solver().area();

    const double &dt = solver().currentTimeStep();

    m_totalTime += dt;

    m_accuDissolutionRate += dt*localDissolutionRate;

    double avgDissolutionRate = m_accuDissolutionRate/m_totalTime;
    const double avgDepositionRate = 1.0;

    m_dMu = avgDissolutionRate/avgDepositionRate;

}
