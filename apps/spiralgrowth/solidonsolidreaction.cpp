#include "solidonsolidreaction.h"

#include "solidonsolidsolver.h"

#include "pressurewall.h"


uint SolidOnSolidReaction::nNeighbors() const
{
    return m_solver.nNeighbors(m_x, m_y);
}

double DiffusionDeposition::calculateDiffusionRate() const
{
    const double &Ew = solver().localPressure(x(), y());
    const double E = (int)nNeighbors() + (int)solver().dim() - 5 + Ew;

    return exp(-solver().alpha()*E - solver().mu());
}

double DiffusionDeposition::calculateDepositionRate() const
{
    if (solver().shadowing())
    {
        return solver().shadowScale(nNeighbors());
    }
    else
    {
        return 1.0;
    }
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

    const uint leftSite = solver().leftSite(x());
    const uint rightSite = solver().rightSite(x());
    const uint bottomSite = solver().bottomSite(y());
    const uint topSite = solver().topSite(y());

    DiffusionDeposition &leftReaction = solver().reaction(leftSite, y());
    DiffusionDeposition &rightReaction = solver().reaction(rightSite, y());
    DiffusionDeposition &bottomReaction = solver().reaction(x(), bottomSite);
    DiffusionDeposition &topReaction = solver().reaction(x(), topSite);

    if (pressureWallEvent.hasStarted())
    {
        pressureWallEvent.registerHeightChange(x(), y());

        pressureWallEvent.findNewHeight();

        pressureWallEvent.recalculateLocalPressure(x(), y());
        pressureWallEvent.recalculateLocalPressure(leftSite, y());
        pressureWallEvent.recalculateLocalPressure(rightSite, y());
        pressureWallEvent.recalculateLocalPressure(x(), bottomSite);
        pressureWallEvent.recalculateLocalPressure(x(), topSite);

        for (uint x = 0; x < solver().length(); ++x)
        {
            for (uint y = 0; y < solver().width(); ++y)
            {

                DiffusionDeposition &reaction = solver().reaction(x, y);

                if (!( (&reaction == this)
                       || (&reaction == &leftReaction)
                       || (&reaction == &rightReaction)
                       || (&reaction == &bottomReaction)
                       || (&reaction == &topReaction)
                       ) )
                {
                    pressureWallEvent.updateRatesFor(reaction);
                }
            }
        }
    }

    calculateRate();
    leftReaction.calculateRate();
    rightReaction.calculateRate();
    bottomReaction.calculateRate();
    topReaction.calculateRate();

}

double DiffusionDeposition::rateExpression()
{
    m_depositionRate = calculateDepositionRate();
    m_diffusionRate = calculateDiffusionRate();

    return m_depositionRate + m_diffusionRate;
}

