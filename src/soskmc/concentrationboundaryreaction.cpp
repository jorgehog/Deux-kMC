#include "concentrationboundaryreaction.h"
#include "sossolver.h"
#include "Events/diffusion/diffusion.h"
#include "Events/confiningsurface/confiningsurface.h"


ConcentrationBoundaryReaction::ConcentrationBoundaryReaction(const uint dim, const uint orientation, SOSSolver &solver) :
    Reaction(),
    m_dim(dim),
    m_orientation(orientation),
    m_location(orientation == 0 ? 0 : (dim == 0 ? solver.length() : solver.width())),
    m_span(dim == 0 ? solver.length() : solver.width()),
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
    return AXTRANS(solver().height, n);
}

void ConcentrationBoundaryReaction::getFreeBoundarSite(const uint n, uint &xi, int &z) const
{
    BADAss(n, <, freeBoundarySites(), "invalid site");

    //unoptimized
    double surfaceHeight = solver().confiningSurfaceEvent().height();
    int surfaceHeightInt = (int)surfaceHeight;

    bool blocked;
    uint nCount = 0;
    for (xi = 0; xi < span(); ++xi)
    {
        for (z = heightAtBoundary(xi) + 1; z < surfaceHeightInt - 1; ++z)
        {
            if (dim() == 0)
            {
                blocked = solver().diffusionEvent().isBlockedPosition(xi, location(), z);
            }

            else
            {
                blocked = solver().diffusionEvent().isBlockedPosition(location(), xi, z);
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
    int surfaceHeightInt = (int)surfaceHeight;

    return (surfaceHeight - surfaceHeightInt);

}

uint ConcentrationBoundaryReaction::freeBoundarySites(bool spam) const
{
    const double surfaceHeight = solver().confiningSurfaceEvent().height();
    const int surfaceHeightInt = (int)surfaceHeight;

    bool blocked;

    uint nSites = 0;

    for (uint xi = 0; xi < span(); ++xi)
    {
        for (int z = heightAtBoundary(xi) + 1; z < surfaceHeightInt - 1; ++z)
        {

            if (dim() == 0)
            {
                blocked = solver().diffusionEvent().isBlockedPosition(xi, location(), z);
            }

            else
            {
                blocked = solver().diffusionEvent().isBlockedPosition(location(), xi, z);
            }

            if (!blocked)
            {
                if (spam)
                {
                    cout << "(" << xi << " " << z << ") ";
                }
                nSites += 1;
            }
        }
    }

    if (spam)
    {
        cout << "\n\n\n" << endl;
    }

    return nSites;
}

bool ConcentrationBoundaryReaction::pointIsOnBoundary(const uint x, const uint y) const
{
    if (dim() == 0)
    {
        return y == location();
    }

    else
    {
        return x == location();
    }
}

double ConcentrationBoundaryReaction::_rateExpression() const
{
    return freeBoundaryArea()/(1-exp(solver().gamma() - 2*solver().alpha()));
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
    return _rateExpression();
}
