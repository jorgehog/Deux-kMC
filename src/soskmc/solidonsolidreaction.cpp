#include "solidonsolidreaction.h"

#include "solidonsolidsolver.h"

#include "Events/pressurewall.h"

#include "Events/cavitydiffusion.h"

#include "../kmcsolver/boundary/boundary.h"


uint SolidOnSolidReaction::nNeighbors() const
{
    return m_solver.nNeighbors(m_x, m_y);
}

double DiffusionDeposition::calculateDiffusionRate() const
{
    const double &Ew = solver().localPressure(x(), y());
    const double E = (int)nNeighbors() + (int)solver().dim() - 5 + Ew;

    return std::exp(-solver().alpha()*E - solver().gamma());
}

double DiffusionDeposition::calculateDepositionRate() const
{
    return solver().localSurfaceSupersaturation(x(), y());
}

bool DiffusionDeposition::isAllowed() const
{
    return true;
}

void DiffusionDeposition::executeAndUpdate()
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

double DiffusionDeposition::rateExpression()
{
    m_depositionRate = calculateDepositionRate();
    m_diffusionRate = calculateDiffusionRate();

    return m_depositionRate + m_diffusionRate;
}
