#include "extraneighbor.h"

ExtraNeighbor::ExtraNeighbor(SOSSolver &solver) :
    LocalPotential(solver),
    m_potentialValues(solver.length(), solver.width(), fill::zeros)
{

}

double ExtraNeighbor::energyFunction(const double dh) const
{
    if (dh > 2)
    {
        return 0;
    }

    else
    {
        return 5*(2-dh);
    }
}

double ExtraNeighbor::energyFunction(const uint x, const uint y) const
{
    const double &h = solver().confiningSurfaceEvent().height();
    const int &hi = solver().height(x, y);

    return energyFunction(h - hi);
}

void ExtraNeighbor::initializeObserver(const Subjects &subject)
{
    if (subject == Subjects::SOLVER)
    {
        for (uint x = 0; x < solver().length(); ++x)
        {
            for (uint y = 0; y < solver().width(); ++y)
            {
                m_potentialValues(x, y) = energyFunction(x, y);
            }
        }
    }
}

void ExtraNeighbor::notifyObserver(const Subjects &subject)
{
    if (subject == Subjects::SOLVER)
    {
        const CurrentSurfaceChange &csc = solver().currentSurfaceChange();

        const uint &x = csc.x;
        const uint &y = csc.y;

        if (csc.type == ChangeTypes::Single)
        {
            const int &value = csc.value;
            SurfaceReaction &r = solver().surfaceReaction(x, y);

            if (value < 1)
            {
                r.setEscapeRate(r.escapeRate() - m_potentialValues(x, y));
                m_potentialValues(x, y) = 0;
            }

            else
            {
                BADAssEqual(m_potentialValues(x, y), 0);
                m_potentialValues(x, y) = energyFunction(x, y);
                r.setEscapeRate(r.escapeRate() + m_potentialValues(x, y));
            }
        }

        else
        {
            BADAssBreak("Extra neighbors and surface diff is not supported.");
        }
    }

    //confining surface
    else
    {
        for (uint x = 0; x < solver().length(); ++x)
        {
            for (uint y = 0; y < solver().width(); ++y)
            {
                SurfaceReaction &r = solver().surfaceReaction(x, y);
                const double prevRate = m_potentialValues(x, y);

                m_potentialValues(x, y) = energyFunction(x, y);

                r.setEscapeRate(r.escapeRate() - prevRate + m_potentialValues(x, y));
            }
        }
    }

#ifndef NDEBUG
    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            BADAssClose(m_potentialValues(x, y), energyFunction(x, y), 1E-10);
        }
    }
#endif

}

double ExtraNeighbor::energy(const uint x, const uint y) const
{
    return m_potentialValues(x, y);
}
