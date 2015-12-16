#include "eqmu.h"

#include "dissolutiondeposition.h"

#include "miscevents.h"


void EqGamma::initialize()
{
    resetCounters();

    update();
}

void EqGamma::execute()
{
    setValue(dGamma());
}

void EqGamma::reset()
{
    update();
}

void EqGamma::restart()
{
    resetCounters();
}

void EqGamma::update()
{
    double localDissolutionRate = 0;

    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            localDissolutionRate += solver().surfaceReaction(x, y).escapeRate();
        }
    }

    localDissolutionRate /= solver().area();

    const double &dt = solver().currentTimeStep();

    m_totalTime += dt;

    m_accuDissolutionRate += dt*localDissolutionRate;

    double avgDissolutionRate = m_accuDissolutionRate/m_totalTime;
    const double avgDepositionRate = 1.0;

    m_dGamma = avgDissolutionRate/avgDepositionRate;

}
