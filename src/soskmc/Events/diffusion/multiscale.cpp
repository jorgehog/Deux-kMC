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

//    const double V = solver().area()*(solver().confiningSurfaceEvent().height() - zMin);
//    const uint N = V*solver().concentration();

//    initializeParticleMatrices(N);

}

void Multiscale::executeFluxBoundaryReaction(const uint x, const uint y, const double z)
{
    (void) x;
    (void) y;
    (void) z;
}

void Multiscale::dump(const uint frameNumber, const string path, const string ext) const
{
    Diffusion::dump(frameNumber, path, ext);
    LatticeDiffusion::dumpDiffusingParticles(frameNumber, path, ext);
    OfflatticeMonteCarlo::dumpDiffusingParticles(frameNumber, path, ext);
}

double Multiscale::concentration() const
{
    return 0;
}

bool Multiscale::hasDiscreteParticles() const
{
    return true;
}

uint Multiscale::numberOfParticles() const
{
    return LatticeDiffusion::numberOfParticles() + OfflatticeMonteCarlo::numberOfParticles();
}

void Multiscale::insertRandomParticle()
{

}

void Multiscale::removeRandomParticle()
{

}

double Multiscale::calculateLocalRateOverD(const uint x, const uint y, const uint n) const
{
    (void) x;
    (void) y;
    (void) n;

    return 0;
}
