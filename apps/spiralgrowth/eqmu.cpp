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
            localDissolutionRate += solver().reaction(x, y).diffusionRate();
        }
    }

    localDissolutionRate /= solver().area();

    const double &dt = solver().currentTimeStep();

    m_totalTime += dt;

    m_accuDissolutionRate += dt*localDissolutionRate;

    double avgDissolutionRate = m_accuDissolutionRate/m_totalTime;

    double inverseKStar;

    if (solver().shadowing())
    {
        m_accuNeighbours += dt*dependency<NNeighbors>("nNeighbors")->getLocalValue();
        const double avgNeighbors = m_accuNeighbours/m_totalTime;

        if (avgNeighbors != avgNeighbors)
        {
            cout << "n " << avgNeighbors << " " << m_accuNeighbours << " " << m_totalTime << " " << dt << endl;
            exit(1);
        }

        inverseKStar = 1./solver().shadowScale(avgNeighbors);
    }
    else
    {
        inverseKStar = 1.;
    }

    m_dMu = avgDissolutionRate*inverseKStar;

    if (m_dMu < 0)
    {
        cout << m_dMu << endl;
        cout << inverseKStar << endl;
        cout << localDissolutionRate << endl;
        cout << avgDissolutionRate << endl;
        exit(1);
    }
}

