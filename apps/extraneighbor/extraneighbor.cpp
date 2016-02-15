#include "extraneighbor.h"

ExtraNeighbor::ExtraNeighbor(SOSSolver &solver, const double energyShift) :
    LocalCachedPotential(solver),
    m_energyShift(energyShift)
{
    solver.registerObserver(this);
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
        const double dh6 = dh2*dh2*dh2;

        return (3*dh + 64/dh6 - 7)*(1 + m_energyShift)/60;
    }
}

void ExtraNeighbor::notifyObserver(const Subjects &subject)
{
    if (subject == Subjects::SOLVER)
    {
        const CurrentSurfaceChange &csc = solver().currentSurfaceChange();

        const uint &x = csc.x;
        const uint &y = csc.y;

        SurfaceReaction &r = solver().surfaceReaction(x, y);

        if (csc.type == ChangeTypes::Single)
        {
            const int &value = csc.value;

            if (value == -1)
            {
                r.setEscapeRate(r.escapeRate() - m_potentialValues(x, y));
                m_potentialValues(x, y) = 0;
            }

            else
            {
                double prevRate = m_potentialValues(x, y);
                m_potentialValues(x, y) = potentialFunction(x, y);

                r.setEscapeRate(r.escapeRate() - prevRate + m_potentialValues(x, y));
            }   
        }

        else
        {
            const uint &xEnd = csc.x1;
            const uint &yEnd = csc.y1;

            SurfaceReaction &rEnd = solver().surfaceReaction(xEnd, yEnd);

            r.setEscapeRate(r.escapeRate() - m_potentialValues(x, y));
            m_potentialValues(x, y) = 0;

            double prevRate = m_potentialValues(xEnd, yEnd);
            m_potentialValues(xEnd, yEnd) = potentialFunction(xEnd, yEnd);
            rEnd.setEscapeRate(rEnd.escapeRate() - prevRate + m_potentialValues(xEnd, yEnd));
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
