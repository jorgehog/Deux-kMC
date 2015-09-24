#include "multiscale.h"

#include "../confiningsurface/confiningsurface.h"

#include "../kmcsolver/boundary/boundary.h"

#include "sosdiffusionreaction.h"

OfflatticeMonteCarloBoundary::OfflatticeMonteCarloBoundary(SOSSolver &solver,
                                                           const double dt,
                                                           const uint boundarySpacing) :
    Diffusion(solver, "MCDiffBoundary"),
    OfflatticeMonteCarlo(solver, dt),
    LatticeDiffusion(solver),
    m_boundarySpacing(boundarySpacing)
{

}

OfflatticeMonteCarloBoundary::~OfflatticeMonteCarloBoundary()
{

}

bool OfflatticeMonteCarloBoundary::checkIfEnoughRoom() const
{
    //derp
    return true;

    const int max = solver().heights().max();
    const double &confinedSurfaceHeight = solver().confiningSurfaceEvent().height();

    return (confinedSurfaceHeight - max >= 2*m_boundarySpacing);
}

void OfflatticeMonteCarloBoundary::execute()
{

}

void OfflatticeMonteCarloBoundary::executeDiffusionReaction(SOSDiffusionReaction *reaction, const int x, const int y, const int z)
{

}


void OfflatticeMonteCarloBoundary::setupInitialConditions()
{
    const double zMin = solver().heights().max() + m_boundarySpacing;

    BADAss(zMin, <, solver().confiningSurfaceEvent().height(), "not enough room for continuum solver.");

    const double V = solver().area()*(solver().confiningSurfaceEvent().height() - zMin);
    const uint N = V*solver().concentration();

    initializeParticleMatrices(N, zMin);

}

void OfflatticeMonteCarloBoundary::executeConcentrationBoundaryReaction(ConcentrationBoundaryReaction *reaction)
{
    //STUFF
}

void OfflatticeMonteCarloBoundary::dump(const uint frameNumber, const string path) const
{
    Diffusion::dump(frameNumber, path);
    LatticeDiffusion::dumpDiffusingParticles(frameNumber, path);
    OfflatticeMonteCarlo::dumpDiffusingParticles(frameNumber, path);
}

double OfflatticeMonteCarloBoundary::calculateLocalRate(const uint x, const uint y, const uint n) const
{
}
