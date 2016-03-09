#include "dissolutiondeposition.h"

#include "sossolver.h"

#include "Events/diffusion/diffusion.h"

#include "Events/confiningsurface/confiningsurface.h"

#include "localpotential.h"


SurfaceReaction::SurfaceReaction(const uint x, const uint y, SOSSolver &solver) :
    SOSReaction(x, y, solver),
    m_isFrozen(false)
{

}

double SurfaceReaction::calculateEscapeRate() const
{
    if (m_isFrozen)
    {
        BADAss(solver().height(x(), y()), >=, m_freezeHeight);

        if (solver().height(x(), y()) == m_freezeHeight)
        {
            return 0;
        }
    }

    if (solver().confiningSurfaceEvent().hasSurface())
    {
        if (solver().height(x(), y()) > solver().confiningSurfaceEvent().height() - 2)
        {
            //no paths
            if (solver().nNeighbors(x(), y()) == 5)
            {
                return 0;
            }
        }
    }

    //we always have the neighbor interaction
    double E = nNeighbors();

    for (const LocalPotential *localPotential : solver().localPotentials())
    {
        E += localPotential->potential(x(), y());
    }

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

    const double escapeRateSingle = std::exp(-solver().alpha()*(E+shift) - gammaTerm);

    if (!solver().diffusionEvent().countPaths())
    {
        return escapeRateSingle;
    }

    else
    {
        if (solver().surfaceDiffusion())
        {
            return solver().numberOfSurroundingSites(x(), y())*escapeRateSingle;
        }

        else
        {
            return solver().numberOfSurroundingSolutionSites(x(), y())*escapeRateSingle;
        }
    }
}

double SurfaceReaction::calculateDepositionRate() const
{
    if (solver().confiningSurfaceEvent().hasSurface())
    {
        if (solver().height(x(), y()) + 1 > solver().confiningSurfaceEvent().height() - 1)
        {
            return 0;
        }
    }

    BADAssBool(solver().diffusionEvent().hasStarted());

    return solver().diffusionEvent().depositionRate(x(), y());
}

void SurfaceReaction::getEscapePath(const uint path, int &dx, int &dy, int &dz) const
{
    bool connectedLeft, connectedRight,
            connectedBottom, connectedTop;

    solver().findConnections(x(), y(),
                             connectedLeft,
                             connectedRight,
                             connectedBottom,
                             connectedTop,
                             false);
    uint n = 0;

    if (!connectedLeft)
    {
        if (n == path)
        {
            dx = -1;
            dy = 0;
            dz = 0;

            return;
        }

        n++;
    }

    if (!connectedRight)
    {
        if (n == path)
        {
            dx = 1;
            dy = 0;
            dz = 0;

            return;
        }

        n++;
    }

    if (!connectedBottom)
    {
        if (n == path)
        {
            dx = 0;
            dy = -1;
            dz = 0;

            return;
        }

        n++;
    }

    if (!connectedTop)
    {
        if (n == path)
        {
            dx = 0;
            dy = 1;
            dz = 0;

            return;
        }

        n++;
    }

    //this leaves up the only option
    dx = 0;
    dy = 0;
    dz = 1;

    BADAssBool(!solver().isBlockedPosition(x(), y(), solver().height(x(), y()) + 1));

}

void SurfaceReaction::freeze()
{
    m_isFrozen = true;
    m_freezeHeight = solver().height(x(), y());
}

void SurfaceReaction::executeAndUpdate()
{
    double r = rate()*rng.uniform();

    if (r < m_depositionRate)
    {
        solver().registerHeightChange(x(), y(), +1);
    }
    else
    {

        if (!solver().surfaceDiffusion())
        {
            solver().registerHeightChange(x(), y(), -1);
        }

        else
        {
            //now we have to decide weather it is dissolution
            //or surface diffusion

            int dx, dy, dz;

            const uint nPaths = solver().numberOfSurroundingSites(x(), y());

            BADAss(nPaths, !=, 0u);

            const uint path = rng.uniform()*nPaths;

            getEscapePath(path, dx, dy, dz);

            int xNew, yNew;
            const int zNew = solver().height(x(), y()) + dz;

            solver().boundaryLatticeTransform(xNew, yNew,
                                              x() + dx,
                                              y() + dy,
                                              zNew);

            if ((dz == 0) && solver().isSurfaceSite(xNew, yNew, zNew))
            {
                solver().registerSurfaceTransition(x(), y(), xNew, yNew);
            }

            else
            {
                solver().registerHeightChange(x(), y(), -1);
            }
        }
    }

}

double SurfaceReaction::rateExpression()
{
    m_depositionRate = calculateDepositionRate();
    m_escapeRate = calculateEscapeRate();

    return m_depositionRate + m_escapeRate;
}

