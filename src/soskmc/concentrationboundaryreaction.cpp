#include "concentrationboundaryreaction.h"
#include "sossolver.h"
#include "Events/diffusion/diffusion.h"
#include "Events/confiningsurface/confiningsurface.h"


ConcentrationBoundaryReaction::ConcentrationBoundaryReaction(const uint dim, const uint orientation, SOSSolver &solver) :
    Reaction(),
    m_dim(dim),
    m_orientation(orientation),
    m_location(orientation == 0 ? 0 : (dim == 0 ? solver.length() : solver.width())),
    m_span(dim == 0 ? solver.width() : solver.length()),
    m_solver(solver)
{
    //location = 0 if orientation = 0 for both dims. If location = 1 we are at the end of the respective dim.
}

ConcentrationBoundaryReaction::~ConcentrationBoundaryReaction()
{

}

double ConcentrationBoundaryReaction::freeBoundaryArea() const
{
    return topFilling()*span() + freeBoundarySites();
}

double ConcentrationBoundaryReaction::dh(const uint n) const
{
    const double &surfaceHeight = solver().confiningSurfaceEvent().height();

    return surfaceHeight - heightAtBoundary(n);
}

const int &ConcentrationBoundaryReaction::heightAtBoundary(const uint n) const
{
    if (dim() == 0)
    {
        return solver().height(location(), n);
    }

    else
    {
        return solver().height(n, location());
    }
}

void ConcentrationBoundaryReaction::getFreeBoundarSite(const uint n, uint &xi, int &z) const
{
    BADAss(n, <, freeBoundarySites(), "invalid site");

    //unoptimized
    double surfaceHeight = solver().confiningSurfaceEvent().height();
    int surfaceHeightInt = floor(surfaceHeight);

    bool blocked;
    uint nCount = 0;
    for (xi = 0; xi < span(); ++xi)
    {
        for (z = heightAtBoundary(xi) + 1; z < surfaceHeightInt; ++z)
        {
            if (dim() == 0)
            {
                blocked = solver().diffusionEvent().isBlockedPosition(location(), xi, z);
            }

            else
            {
                blocked = solver().diffusionEvent().isBlockedPosition(xi, location(), z);
            }

            if (!blocked)
            {
                if (nCount == n)
                {
                    return;
                }

                nCount++;
            }
        }
    }
}

double ConcentrationBoundaryReaction::topFilling() const
{
    double surfaceHeight = solver().confiningSurfaceEvent().height();
    int surfaceHeightInt = floor(surfaceHeight);

    return (surfaceHeight - surfaceHeightInt);

}

uint ConcentrationBoundaryReaction::freeBoundarySites() const
{
    const double surfaceHeight = solver().confiningSurfaceEvent().height();
    const int surfaceHeightInt = (const int)surfaceHeight;

    bool blocked;

    uint nSites = 0;

    for (uint xi = 0; xi < span(); ++xi)
    {
        for (int z = heightAtBoundary(xi) + 1; z < surfaceHeightInt; ++z)
        {

            if (dim() == 0)
            {
                blocked = solver().diffusionEvent().isBlockedPosition(location(), xi, z);
            }

            else
            {
                blocked = solver().diffusionEvent().isBlockedPosition(xi, location(), z);
            }

            if (!blocked)
            {
                nSites += 1;
            }
        }
    }

    return nSites;
}

bool ConcentrationBoundaryReaction::pointIsOnBoundary(const uint x, const uint y) const
{
    if (dim() == 0)
    {
        return x == location();
    }

    else
    {
        return y == location();
    }
}

double ConcentrationBoundaryReaction::_rateExpression(const double freeArea) const
{
    return freeArea/(1-solver().concentration());
}

bool ConcentrationBoundaryReaction::isAllowed() const
{
    return freeBoundaryArea() != 0;
}

void ConcentrationBoundaryReaction::executeAndUpdate()
{
    solver().diffusionEvent().executeConcentrationBoundaryReaction(this);

    calculateRate();
}

double ConcentrationBoundaryReaction::rateExpression()
{
    return _rateExpression(freeBoundaryArea());
}


void ConcentrationBoundaryReaction::initializeObserver(const Subjects &subject)
{
    (void) subject;
}

void ConcentrationBoundaryReaction::notifyObserver(const Subjects &subject)
{
    if (subject == Subjects::SOLVER)
    {
        const uint &x = solver().currentSurfaceChange().x;
        const uint &y = solver().currentSurfaceChange().y;

        bool onBoundary = pointIsOnBoundary(x, y);

        if (solver().currentSurfaceChange().type == ChangeTypes::Double)
        {
            const uint &x1 = solver().currentSurfaceChange().x1;
            const uint &y1 = solver().currentSurfaceChange().y1;

            onBoundary = onBoundary || pointIsOnBoundary(x1, y1);
        }

        if (onBoundary)
        {
            m_solver.registerAffectedReaction(this);
        }
    }

    else
    {
        const double heightChange = solver().confiningSurfaceEvent().height() -
                solver().confiningSurfaceEvent().currentConfinementChange().prevHeight;
        changeRate(rate() + _rateExpression(heightChange*span()));
    }
}
