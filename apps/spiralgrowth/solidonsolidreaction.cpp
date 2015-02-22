#include "solidonsolidreaction.h"

#include "solidonsolidsolver.h"


uint SolidOnSolidReaction::nNeighbors() const
{
    return m_solver.nNeighbors(m_x, m_y);
}


bool DiffusionDeposition::isAllowed() const
{
    return true;
}

void DiffusionDeposition::executeAndUpdate()
{
    double r = rate()*rng.uniform();

    if (r < 1)
    {
        system().registerHeightChange(x(), y(), +1);
    }
    else
    {
        system().registerHeightChange(x(), y(), -1);
    }

    calculateRate();
    system().reaction(system().leftSite(x()), y()).calculateRate();
    system().reaction(system().rightSite(x()), y()).calculateRate();
    system().reaction(x(), system().bottomSite(y())).calculateRate();
    system().reaction(x(), system().topSite(y())).calculateRate();

}

double DiffusionDeposition::rateExpression()
{
    return 1 + exp(-system().alpha()*((int)nNeighbors() - 4));
}

