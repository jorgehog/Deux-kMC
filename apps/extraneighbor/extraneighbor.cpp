#include "extraneighbor.h"

ExtraNeighbor::ExtraNeighbor(SOSSolver &solver) :
    LocalCachedPotential(solver)
{

}

double ExtraNeighbor::energyFunction(const double dh) const
{
    BADAss(dh, >=, 1);

    if (dh > 2)
    {
        return 0;
    }

    else
    {
        const double dh2 = dh*dh;
        return m_scaling*(1/(dh2*dh2*dh2) - m_shift);
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
                m_potentialValues(x, y) = potentialFunction(x, y);
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

                m_potentialValues(x, y) = potentialFunction(x, y);

                r.setEscapeRate(r.escapeRate() - prevRate + m_potentialValues(x, y));
            }
        }
    }

}

double ExtraNeighbor::potentialFunction(const uint x, const uint y) const
{
    const double &h = solver().confiningSurfaceEvent().height();
    const int &hi = solver().height(x, y);

    return energyFunction(h - hi);
}
