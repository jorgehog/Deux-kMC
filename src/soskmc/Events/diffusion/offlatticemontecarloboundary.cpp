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
    const uint THRESH = 1;
    if (cycle() % THRESH == 0)
    {
        dump(cycle()/THRESH);
    }

    const double &h = solver().confiningSurfaceEvent().height();
    const int zMin = solver().heights().min();

    lammpswriter writer(4, "kmcdiff", "/tmp");
    writer.setSystemSize(solver().length(), solver().width(), h, 0, 0, zMin);
    writer.initializeNewFile(cycle()/THRESH);

    for (uint n = 0; n < nOfflatticeParticles(); ++n)
    {
        uint x = 0;
        uint y = 0;
        int z = 0;

        BADAssBool(!solver().isBlockedPosition(x, y, z));

        writer << 0
               << x
               << y
               << z;
    }


}

void OfflatticeMonteCarloBoundary::initialize()
{
}

void OfflatticeMonteCarloBoundary::reset()
{
}

void OfflatticeMonteCarloBoundary::setupInitialConditions()
{
    const double zMin = solver().heights().max() + m_boundarySpacing;
    const double V = solver().area()*(solver().confiningSurfaceEvent().height() - zMin);
    const uint N = V*solver().concentration();

    initializeParticleMatrices(N, zMin);

    const double latticeVolume = solver().volume() - V;

    const uint nLatticeParticles = latticeVolume*solver().concentration();

    const int hMin = solver().heights().min();

    uint x0;
    uint y0;
    uint z0;

    uint n = 0;
    while (n < nLatticeParticles)
    {
        do
        {
            x0 = rng.uniform()*solver().length();
            y0 = rng.uniform()*solver().width();
            z0 = hMin + rng.uniform()*(zMin - hMin);

        } while(solver().isBlockedPosition(x0, y0, z0));

        cout << "placed at " << x0 << " " << y0 << " " << z0 << endl;

        n++;
    }

}

double OfflatticeMonteCarloBoundary::depositionRate(const uint x, const uint y) const
{
    (void) x;
    (void) y;

    //Deposition is modelled as a diffusion reaction and is not explicitly treated.
    return 0;
}

void OfflatticeMonteCarloBoundary::registerHeightChange(const uint x, const uint y, const int delta)
{
    (void) (x+y+delta);
    BADAssBool(checkIfEnoughRoom());

    //position particle on random surrounding site
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
