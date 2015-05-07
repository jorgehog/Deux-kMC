#include "solidonsolidreaction.h"

#include "solidonsolidsolver.h"

#include "Events/pressurewall.h"

#include "../kmcsolver/boundary/boundary.h"


uint SolidOnSolidReaction::nNeighbors() const
{
    return m_solver.nNeighbors(m_x, m_y);
}

double DiffusionDeposition::calculateDiffusionRate() const
{
    const double &Ew = solver().localPressure(x(), y());
    const double E = (int)nNeighbors() + (int)solver().dim() - 5 + Ew;

    return std::exp(-solver().alpha()*E - solver().mu());
}

double DiffusionDeposition::calculateDepositionRate() const
{
//    if (solver().shadowing())
//    {
//        return solver().shadowScale(nNeighbors());
//    }
//    else
//    {
//        return 1.0;
//    }

    return solver().localSurfaceSupersaturation(x(), y());
}

bool DiffusionDeposition::isAllowed() const
{
    return true;
}

void DiffusionDeposition::executeAndUpdate()
{
    PressureWall &pressureWallEvent = solver().pressureWallEvent();

    double r = rate()*rng.uniform();

    if (r < m_depositionRate)
    {
        solver().registerHeightChange(x(), y(), +1);
    }
    else
    {
        solver().registerHeightChange(x(), y(), -1);
    }

    const uint left = solver().leftSite(x());
    const uint right = solver().rightSite(x());
    const uint bottom = solver().bottomSite(y());
    const uint top = solver().topSite(y());

    uint n = 1;
    m_affectedReactions.at(0) = this;

    if (!solver().boundary(0)->isBlocked(left))
    {
        m_affectedReactions.at(n) = &solver().reaction(left, y());
        n++;
    }

    if (!solver().boundary(0)->isBlocked(right))
    {
        m_affectedReactions.at(n) = &solver().reaction(right, y());
        n++;
    }

    if (!solver().boundary(1)->isBlocked(top))
    {
        m_affectedReactions.at(n) = &solver().reaction(x(), top);
        n++;
    }

    if (!solver().boundary(1)->isBlocked(bottom))
    {
        m_affectedReactions.at(n) = &solver().reaction(x(), bottom);
        n++;
    }

    auto start = m_affectedReactions.begin();
    auto end = start + n;

    if (pressureWallEvent.hasStarted())
    {
        pressureWallEvent.registerHeightChange(x(), y());

        pressureWallEvent.findNewHeight();

        for (uint i = 0; i < n; ++i)
        {
            DiffusionDeposition *reaction = m_affectedReactions.at(i);
            pressureWallEvent.recalculateLocalPressure(reaction->x(), reaction->y());
        }

        for (uint x = 0; x < solver().length(); ++x)
        {
            for (uint y = 0; y < solver().width(); ++y)
            {

                DiffusionDeposition &reaction = solver().reaction(x, y);

                if (std::find(start, end, &reaction) == end)
                {
                    pressureWallEvent.updateRatesFor(reaction);
                }
            }
        }
    }

    for (uint i = 0; i < n; ++i)
    {
        DiffusionDeposition *reaction = m_affectedReactions.at(i);
        reaction->calculateRate();
    }

}

double DiffusionDeposition::rateExpression()
{
    m_depositionRate = calculateDepositionRate();
    m_diffusionRate = calculateDiffusionRate();

    return m_depositionRate + m_diffusionRate;
}

std::vector<DiffusionDeposition*> DiffusionDeposition::m_affectedReactions(5);
