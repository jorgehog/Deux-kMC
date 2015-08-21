#include "concentrationboundaryreaction.h"
#include "sossolver.h"
#include "Events/diffusion/diffusion.h"
#include "Events/confiningsurface/confiningsurface.h"


concentrationBoundaryReaction::concentrationBoundaryReaction(const uint dim, const uint orientation, SOSSolver &solver) :
    Reaction(),
    m_dim(dim),
    m_orientation(orientation),
    m_location(orientation == 0 ? 0 : (dim == 0 ? solver.length() : solver.width())),
    m_span(dim == 0 ? solver.length() : solver.width()),
    m_solver(solver)
{
    //location = 0 if orientation = 0 for both dims. If location = 1 we are at the end of the respective dim.
}

concentrationBoundaryReaction::~concentrationBoundaryReaction()
{

}

double concentrationBoundaryReaction::freeBoundaryArea() const
{
    double surfaceHeight = solver().confiningSurfaceEvent().height();
    int surfaceHeightInt = (int)surfaceHeight;

    double area = span()*(surfaceHeight - surfaceHeightInt);

    for (uint xi = 0; xi < span(); ++xi)
    {
        for (int z = base(xi); z < surfaceHeightInt; ++z)
        {
            if (!AXTRANS(solver().diffusionEvent().isBlockedPosition, xi, z))
            {
                area += 1;
            }
        }
    }

    return area;
}

double concentrationBoundaryReaction::dh(const uint n) const
{
    const double &surfaceHeight = solver().confiningSurfaceEvent().height();

    return surfaceHeight - base(n);
}

const int &concentrationBoundaryReaction::base(const uint n) const
{
    return AXTRANS(solver().height, n);
}

bool concentrationBoundaryReaction::isAllowed() const
{
    return freeBoundaryArea() != 0;
}

void concentrationBoundaryReaction::executeAndUpdate()
{
    //choose z randomly across boundary. Subtract Nb and find correct z.

    uint x = 0;
    uint y = 0;
    int z = 0;
    solver().diffusionEvent().insertDiffusingParticle(x, y, z);
}

double concentrationBoundaryReaction::rateExpression()
{
    return freeBoundaryArea()/(1-exp(solver().gamma() - 2*solver().alpha()));
}
