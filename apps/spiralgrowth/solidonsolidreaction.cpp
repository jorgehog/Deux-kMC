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
        solver().registerHeightChange(x(), y(), +1);
    }
    else
    {
        solver().registerHeightChange(x(), y(), -1);
    }

    calculateRate();
    solver().reaction(solver().leftSite(x()), y()).calculateRate();
    solver().reaction(solver().rightSite(x()), y()).calculateRate();
    solver().reaction(x(), solver().bottomSite(y())).calculateRate();
    solver().reaction(x(), solver().topSite(y())).calculateRate();

}

double DiffusionDeposition::rateExpression()
{
    return 1 + exp(-solver().alpha()*((int)nNeighbors() - 5 + (int)solver().dim()));
}

