#include "spiralgrowthsolver.h"

SolidOnSolidSolver::SolidOnSolidSolver(const uint length, const double alpha) :
    KMCSolver(),
    m_length(length),
    m_alpha(alpha),
    m_heights(length, fill::zeros),
    m_siteReactions(length)
{
    for (uint x = 0; x < length; ++x)
    {
        m_siteReactions(x) = new DiffusionDeposition(x, *this);
    }
}

SolidOnSolidSolver::~SolidOnSolidSolver()
{

}

uint SolidOnSolidSolver::nNeighbors(const uint site) const
{

    bool hasLeftNeighbor = m_heights(leftSite(site)) >= m_heights(site);
    bool hasRightNeighbor = m_heights(rightSite(site)) >= m_heights(site);

    if (hasLeftNeighbor && hasRightNeighbor)
    {
        return 3;
    }

    else if (hasLeftNeighbor || hasRightNeighbor)
    {
        return 2;
    }

    else
    {
        return 1;
    }
}

uint SolidOnSolidSolver::leftSite(const uint site, const uint n) const
{
    return (site + length() - n)%length();
}

uint SolidOnSolidSolver::rightSite(const uint site, const uint n) const
{
    return (site + n)%length();
}

bool DiffusionDeposition::isAllowed() const
{
    return true;
}

void DiffusionDeposition::executeAndUpdate()
{
    double r = (1 + rate())*rng.uniform();

    if (r < 1)
    {
        system().registerHeightChange(x(), +1);
    }
    else
    {
        system().registerHeightChange(x(), -1);
    }

    calculateRate();
    system().reaction(system().leftSite(x())).calculateRate();
    system().reaction(system().rightSite(x())).calculateRate();
}

double DiffusionDeposition::rateExpression()
{
    return 1 + exp(-system().alpha()*(nNeighbors() - 2));
}


uint SolidOnSolidSolver::numberOfReactions() const
{
    return m_siteReactions.size();
}

Reaction *SolidOnSolidSolver::getReaction(const uint n) const
{
    return m_siteReactions(n);
}
