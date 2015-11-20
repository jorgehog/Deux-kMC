#include "multiscale.h"

#include "../confiningsurface/confiningsurface.h"

#include "../kmcsolver/boundary/boundary.h"

#include "sosdiffusionreaction.h"

#include "sossolver.h"

Multiscale::Multiscale(SOSSolver &solver,
                       const double dt,
                       const uint boundarySpacing) :
    Diffusion(solver, "MCDiffBoundary"),
    OfflatticeMonteCarlo(solver, dt),
    LatticeDiffusion(solver),
    m_boundarySpacing(boundarySpacing)
{

}

Multiscale::~Multiscale()
{

}

bool Multiscale::checkIfEnoughRoom() const
{
    //derp
    return true;

    const int max = solver().heights().max();
    const double &confinedSurfaceHeight = solver().confiningSurfaceEvent().height();

    return (confinedSurfaceHeight - max >= 2*m_boundarySpacing);
}

void Multiscale::execute()
{

}

void Multiscale::executeDiffusionReaction(SOSDiffusionReaction *reaction, const int x, const int y, const int z)
{
    (void) reaction;
    (void) x;
    (void) y;
    (void) z;
}


void Multiscale::initializeObserver(const Subjects &subject)
{
    (void) subject;

    const double zMin = solver().heights().max() + m_boundarySpacing;

    BADAss(zMin, <, solver().confiningSurfaceEvent().height(), "not enough room for continuum solver.");

    const double V = solver().area()*(solver().confiningSurfaceEvent().height() - zMin);
    const uint N = V*solver().concentration();

    initializeParticleMatrices(N, zMin);

}

void Multiscale::executeConcentrationBoundaryReaction(ConcentrationBoundaryReaction *reaction)
{
    (void) reaction;
}

void Multiscale::dump(const uint frameNumber, const string path) const
{
    Diffusion::dump(frameNumber, path);
    LatticeDiffusion::dumpDiffusingParticles(frameNumber, path);
    OfflatticeMonteCarlo::dumpDiffusingParticles(frameNumber, path);
}

double Multiscale::concentration() const
{
    return 0;
}

double Multiscale::calculateLocalRateOverD(const uint x, const uint y, const uint n) const
{
    (void) x;
    (void) y;
    (void) n;

    return 0;
}
