
#include "dissolutiondeposition.h"

#include "sossolver.h"

#include "Events/diffusion/diffusion.h"

#include "Events/confiningsurface/confiningsurface.h"

double DissolutionDeposition::calculateDissolutionRate() const
{
    if (solver().confiningSurfaceEvent().hasSurface())
    {
        if (solver().height(x(), y()) >= solver().confiningSurfaceEvent().height() - 1)
        {
            uint nnBelow = solver().calculateNNeighbors(x(), y(), solver().height(x(), y()) - 1);

            //no paths
            if (nnBelow == 5)
            {
                return 0;
            }
        }
    }

    const double &Ew = solver().confinementEnergy(x(), y());
    const double E = (int)nNeighbors() + (int)solver().surfaceDim() - 5 + Ew;

    double gammaTerm;

    if (solver().concentrationIsZero())
    {
        gammaTerm = 0;
    }

    else
    {
        gammaTerm = solver().gamma();
    }

    return solver().diffusionEvent().dissolutionPaths(x(), y())*std::exp(-solver().alpha()*E - gammaTerm);
}

double DissolutionDeposition::calculateDepositionRate() const
{
    //derp when wall moves away, this reaction should have its rate changed back. Very rare?
    if (solver().confiningSurfaceEvent().hasSurface())
    {
        if (solver().height(x(), y()) + 1 > solver().confiningSurfaceEvent().height() - 1)
        {
            return 0;
        }
    }

    if (!solver().diffusionEvent().hasStarted())
    {
        return 1.0;
    }

    if (solver().concentrationIsZero())
    {
        return 0;
    }

    return solver().diffusionEvent().depositionRate(x(), y());
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

