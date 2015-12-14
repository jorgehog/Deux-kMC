
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

    //this ensures that the expression goes as exp(-alpha*(n - 2 or 3))
    const double E = (int)nNeighbors() + Ew;

    double gammaTerm;
    double shift;


    if (solver().concentrationIsZero())
    {
        //we use an expression equal to exp(-alpha*n)/c0
        gammaTerm = 0; //now he have divided by nothing
        shift = 0;
    }

    else
    {
        //now we have divided by the concentration/c0
        gammaTerm = solver().gamma();
        shift = (int)solver().surfaceDim() - 5;
    }

    return solver().diffusionEvent().dissolutionPaths(x(), y())*std::exp(-solver().alpha()*(E+shift) - gammaTerm);
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
    m_dissolutionRate = calculateDissolutionRate();

    return m_depositionRate + m_dissolutionRate;
}

