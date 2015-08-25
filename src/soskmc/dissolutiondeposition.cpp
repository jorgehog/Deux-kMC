#include "dissolutiondeposition.h"

#include "sossolver.h"
#include "Events/diffusion/diffusion.h"


double DissolutionDeposition::calculateDissolutionRate() const
{
    //rate is zero if blocked by diffusing particles.
    if (solver().numberOfSurroundingSolutionSites(x(), y()) == 0)
    {
        return 0.0;
    }

    const double &Ew = solver().confinementEnergy(x(), y());
    const double E = (int)nNeighbors() + (int)solver().surfaceDim() - 5 + Ew;

    return solver().numberOfSurroundingSolutionSites(x(), y())*std::exp(-solver().alpha()*E - solver().gamma());
}

double DissolutionDeposition::calculateDepositionRate() const
{
    return solver().depositionRate(x(), y());
}

bool DissolutionDeposition::isAllowed() const
{
    return true;
}

void DissolutionDeposition::executeAndUpdate()
{
    double r = rate()*rng.uniform();

    if (r < m_depositionRate)
    {
        solver().registerHeightChange(x(), y(), +1);
    }
    else
    {
        solver().registerHeightChange(x(), y(), -1);
    }

}

double DissolutionDeposition::rateExpression()
{
    m_depositionRate = calculateDepositionRate();
    m_diffusionRate = calculateDissolutionRate();

    return m_depositionRate + m_diffusionRate;
}

