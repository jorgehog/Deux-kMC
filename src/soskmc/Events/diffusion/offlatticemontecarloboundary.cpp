#include "offlatticemontecarloboundary.h"

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
    dump(cycle());
}


void OfflatticeMonteCarloBoundary::setupInitialConditions()
{
    const double zMin = solver().heights().max() + m_boundarySpacing;

    BADAss(zMin, <, solver().confiningSurfaceEvent().height(), "not enough room for continuum solver.");

    const double V = solver().area()*(solver().confiningSurfaceEvent().height() - zMin);
    const uint N = V*solver().concentration();

    initializeParticleMatrices(N, zMin);

    //derp
//    const double latticeVolume = solver().volume() - V;

//    const uint nLatticeParticles = latticeVolume*solver().concentration();

//    const int hMin = solver().heights().min();

//    uint x0;
//    uint y0;
//    uint z0;

//    uint n = 0;
//    while (n < nLatticeParticles)
//    {
//        do
//        {
//            x0 = rng.uniform()*solver().length();
//            y0 = rng.uniform()*solver().width();
//            z0 = hMin + rng.uniform()*(zMin - hMin);

//        } while(solver().isBlockedPosition(x0, y0, z0) ||
//                solver().isSurfaceSite(x0, y0, z0) ||
//                isBlockedPosition(x0, y0, z0));

//        addDiffusionReactant(x0, y0, z0, false);

//        n++;
//    }

}

void OfflatticeMonteCarloBoundary::dump(const uint frameNumber) const
{
    Diffusion::dump(frameNumber);
    LatticeDiffusion::dump(frameNumber);
    OfflatticeMonteCarlo::dump(frameNumber);
}

void OfflatticeMonteCarloBoundary::onInsertParticle(const double x, const double y, const double z)
{
    (void) (x + y + z);

    //update flux if on boundary
}

void OfflatticeMonteCarloBoundary::onRemoveParticle(const uint n)
{
    (void) n;

    //update flux if on boundary
}


