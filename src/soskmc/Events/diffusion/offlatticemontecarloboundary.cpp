#include "offlatticemontecarloboundary.h"

#include "../confiningsurface/confiningsurface.h"

#include "../kmcsolver/boundary/boundary.h"

OfflatticeMonteCarloBoundary::OfflatticeMonteCarloBoundary(SOSSolver &solver,
                                                           const double dt,
                                                           const uint boundarySpacing) :
    OfflatticeMonteCarlo(solver, dt, "MCDiffBound"),
    m_boundarySpacing(boundarySpacing)
{

}

OfflatticeMonteCarloBoundary::~OfflatticeMonteCarloBoundary()
{

}

bool OfflatticeMonteCarloBoundary::checkIfEnoughRoom() const
{
    const int max = solver().heights().max();
    const double &confinedSurfaceHeight = solver().confiningSurfaceEvent().height();

    return (confinedSurfaceHeight - max >= 2*m_boundarySpacing);
}

void OfflatticeMonteCarloBoundary::execute()
{
}

void OfflatticeMonteCarloBoundary::initialize()
{
}

void OfflatticeMonteCarloBoundary::reset()
{
}

void OfflatticeMonteCarloBoundary::setupInitialConditions()
{
}

double OfflatticeMonteCarloBoundary::depositionRate(const uint x, const uint y) const
{
    //find out which are on surrounding sites and divide by number of surrounding sites
    return 1.0 + 0*x*y;
}

void OfflatticeMonteCarloBoundary::registerHeightChange(const uint x, const uint y, const int delta)
{
    (void) (x+y+delta);
    BADAssBool(checkIfEnoughRoom());

    //position particle on random surrounding site
}


void OfflatticeMonteCarloBoundary::onInsertParticle(const double x, const double y, const double z)
{
}

void OfflatticeMonteCarloBoundary::onRemoveParticle(const uint n)
{
}
