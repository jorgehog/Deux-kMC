#include "offlatticemontecarloboundary.h"

#include "../confiningsurface/confiningsurface.h"

#include "../kmcsolver/boundary/boundary.h"

#include "sosdiffusionreaction.h"

OfflatticeMonteCarloBoundary::OfflatticeMonteCarloBoundary(SOSSolver &solver,
                                                           const double dt,
                                                           const uint boundarySpacing) :
    OfflatticeMonteCarlo(solver, dt, "MCDiffBound"),
    m_boundarySpacing(boundarySpacing),
    m_mutexSolver(solver)
{

}

OfflatticeMonteCarloBoundary::~OfflatticeMonteCarloBoundary()
{
    deleteQueuedReactions();

    for (SOSDiffusionReaction *reaction : m_diffusionReactions)
    {
        delete reaction;
    }

    m_diffusionReactions.clear();
}

bool OfflatticeMonteCarloBoundary::checkIfEnoughRoom() const
{
    //derp
    return true;

    const int max = solver().heights().max();
    const double &confinedSurfaceHeight = solver().confiningSurfaceEvent().height();

    return (confinedSurfaceHeight - max >= 2*m_boundarySpacing);
}

void OfflatticeMonteCarloBoundary::removeDiffusionReactant(SOSDiffusionReaction *reaction, bool _delete)
{
    m_mutexSolver.removeReaction(reaction);

    auto &r = m_diffusionReactions;
    r.erase( std::remove( r.begin(), r.end(), reaction ), r.end() );

    if (_delete)
    {
        delete reaction; //counters new in addDiffusionReactant
    }

    else
    {
        //dangerous to delete a reaction if it is used somewhere else.
        //Should use smart pointers really.
        m_deleteQueue.push_back(reaction);
    }
}

void OfflatticeMonteCarloBoundary::removeDiffusionReactant(const uint x, const uint y, const int z, bool _delete)
{
    removeDiffusionReactant(diffusionReaction(x, y, z), _delete);
}

SOSDiffusionReaction *OfflatticeMonteCarloBoundary::diffusionReaction(const uint x, const uint y, const int z) const
{
    auto &r = m_diffusionReactions;
    const auto &res = std::find_if(r.begin(), r.end(), [&x, &y, &z] (const SOSDiffusionReaction *reaction)
    {
        return (x == reaction->x()) && (y == reaction->y()) && (z == reaction->z());
    });

    if (res == r.end())
    {
        return NULL;
    }

    else
    {
        return *res;
    }
}

SOSDiffusionReaction *OfflatticeMonteCarloBoundary::diffusionReaction(const uint n) const
{
    return m_diffusionReactions.at(n);
}

void OfflatticeMonteCarloBoundary::clearDiffusionReactions()
{
    for (SOSDiffusionReaction *reaction : m_diffusionReactions)
    {
        m_mutexSolver.removeReaction(reaction);
        delete reaction;
    }

    m_diffusionReactions.clear();
}

void OfflatticeMonteCarloBoundary::deleteQueuedReactions()
{
    for (SOSDiffusionReaction *reaction : m_deleteQueue)
    {
        delete reaction;
    }

    m_deleteQueue.clear();
}

SOSDiffusionReaction *OfflatticeMonteCarloBoundary::addDiffusionReactant(const uint x, const uint y, const int z, bool setRate)
{
    BADAssBool(!isBlockedPosition(x, y, z), "spot already taken");

    SOSDiffusionReaction *reaction = new SOSDiffusionReaction(m_mutexSolver, x, y, z);

    m_diffusionReactions.push_back(reaction);
    m_mutexSolver.addReaction(reaction);

    if (setRate)
    {
        reaction->calculateRate();
    }

    return reaction;
}

void OfflatticeMonteCarloBoundary::execute()
{
    deleteQueuedReactions();

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

    for(const SOSDiffusionReaction *reaction : m_diffusionReactions)
    {
        const uint &x = reaction->x();
        const uint &y = reaction->y();
        const int &z = reaction->z();

        BADAssBool(!solver().isBlockedPosition(x, y, z), "Illigal particle position", [&] ()
        {
            const int h = solver().height(x, y);
            const double hconf = solver().confiningSurfaceEvent().height();
            BADAssSimpleDump(x, y, z, h, hconf);
        });

        writer << 3
               << x
               << y
               << z;
    }

    writer.finalize();

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

    BADAss(zMin, <, solver().confiningSurfaceEvent().height(), "not enough room for continuum solver.");

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

        } while(solver().isBlockedPosition(x0, y0, z0) ||
                solver().isSurfaceSite(x0, y0, z0) ||
                isBlockedPosition(x0, y0, z0));

        addDiffusionReactant(x0, y0, z0, false);

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

void OfflatticeMonteCarloBoundary::executeDiffusionReaction(SOSDiffusionReaction *reaction,
                                                         const uint x, const uint y, const int z)
{
    int roof = 100000;
    if (z > roof)
    {
        //Derp boundary crossed... fix this
        removeDiffusionReactant(reaction);
        return;
    }

    if (solver().isSurfaceSite(x, y, z))
    {
        m_mutexSolver.registerHeightChange(x, y, 1);
        removeDiffusionReactant(reaction);

        int zAbove = z+1;
        while (isBlockedPosition(x, y, zAbove))
        {
            m_mutexSolver.registerHeightChange(x, y, 1);
            removeDiffusionReactant(x, y, zAbove);
            zAbove++;
        }
    }

    else
    {
        reaction->setX(x);
        reaction->setY(y);
        reaction->setZ(z);
    }
}

bool OfflatticeMonteCarloBoundary::isBlockedPosition(const uint x, const uint y, const int z) const
{
    return diffusionReaction(x, y, z) != NULL;;
}

void OfflatticeMonteCarloBoundary::registerHeightChange(const uint x, const uint y, const int delta)
{
    BADAssBool(checkIfEnoughRoom());

    //position particle on random surrounding site
    //add one to solution site heights because at this stage height(x,y) is already updated
    if (delta == -1)
    {
        const uint randomSite = rng.uniform()*solver().numberOfSurroundingSolutionSites(x, y, solver().height(x, y)+1);

        int dx, dy, dz;
        solver().getSolutionSite(x, y, solver().height(x, y)+1, dx, dy, dz, randomSite);

        const uint xNew = solver().boundary(0)->transformCoordinate((int)x + dx);
        const uint yNew = solver().boundary(1)->transformCoordinate((int)y + dy);
        const int zNew = solver().height(x, y) + dz + 1;

        BADAssBool(!solver().isSurfaceSite(xNew, yNew, zNew), "adding diff reaction to surface.", [&] ()
        {
            int h = solver().height(x, y);
            int hNew = solver().height(xNew, yNew);
            BADAssSimpleDump(x, y, h, xNew, yNew, zNew, hNew);
        });

        addDiffusionReactant(xNew, yNew, zNew);
    }

    //delta == 1 case of deposition is treated from the deposition reactions execute function.
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

