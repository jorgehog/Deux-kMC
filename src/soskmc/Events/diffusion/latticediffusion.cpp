#include "latticediffusion.h"

#include "../confiningsurface/confiningsurface.h"

#include "../kmcsolver/boundary/boundary.h"

#include "sosdiffusionreaction.h"


LatticeDiffusion::LatticeDiffusion(SOSSolver &solver) :
    Diffusion(solver, "latticeDiffusion"),
    m_mutexSolver(solver)
{

}

LatticeDiffusion::~LatticeDiffusion()
{
    deleteQueuedReactions();

    for (SOSDiffusionReaction *reaction : m_diffusionReactions)
    {
        delete reaction;
    }

    m_diffusionReactions.clear();
}


void LatticeDiffusion::removeDiffusionReactant(SOSDiffusionReaction *reaction, bool _delete)
{

    m_mutexSolver.removeReaction(reaction);

#ifndef NDEBUG
    auto comp = [&reaction] (const SOSDiffusionReaction *r)
    {
        return (reaction->x() == r->x()) && (reaction->y() == r->y()) && (reaction->z() == r->z());
    };
#endif

    auto &r = m_diffusionReactions;
    BADAss(std::find_if(r.begin(), r.end(), comp), !=, r.end());
    r.erase( std::remove( r.begin(), r.end(), reaction ), r.end() );
    BADAss(std::find_if(r.begin(), r.end(), comp), ==, r.end());
//    cout << "removed " << reaction->x() << " " << reaction->y() << " " << reaction->z() << endl;
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

void LatticeDiffusion::removeDiffusionReactant(const uint x, const uint y, const int z, bool _delete)
{
    removeDiffusionReactant(diffusionReaction(x, y, z), _delete);
}

SOSDiffusionReaction *LatticeDiffusion::diffusionReaction(const uint x, const uint y, const int z) const
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

SOSDiffusionReaction *LatticeDiffusion::diffusionReaction(const uint n) const
{
    return m_diffusionReactions.at(n);
}

void LatticeDiffusion::clearDiffusionReactions()
{
    for (SOSDiffusionReaction *reaction : m_diffusionReactions)
    {
        m_mutexSolver.removeReaction(reaction);
        delete reaction;
    }

    m_diffusionReactions.clear();
}

void LatticeDiffusion::deleteQueuedReactions()
{
    for (SOSDiffusionReaction *reaction : m_deleteQueue)
    {
        delete reaction;
    }

    m_deleteQueue.clear();
}

void LatticeDiffusion::dump(const uint frameNumber) const
{
    const double &h = solver().confiningSurfaceEvent().height();
    const int zMin = solver().heights().min();

    lammpswriter writer(4, "kmcdiff", "/tmp");
    writer.setSystemSize(solver().length(), solver().width(), h, 0, 0, zMin);
    writer.initializeNewFile(frameNumber);

    for(const SOSDiffusionReaction *reaction : m_diffusionReactions)
    {
        const uint &x = reaction->x();
        const uint &y = reaction->y();
        const int &z = reaction->z();

        writer << 3
               << x
               << y
               << z;
    }

    writer.finalize();
}

SOSDiffusionReaction *LatticeDiffusion::addDiffusionReactant(const uint x, const uint y, const int z, bool setRate)
{

    if (solver().isOutsideBox(x, y))
    {
        //derp fix concentration stuff
        return NULL;
    }

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

void LatticeDiffusion::execute()
{
    deleteQueuedReactions();

    Diffusion::dump(cycle());
    dump(cycle());

#ifndef NDEBUG
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

        BADAssBool(!solver().isSurfaceSite(x, y, z), "Illigal particle position", [&] ()
        {
            const int h = solver().height(x, y);
            const double hconf = solver().confiningSurfaceEvent().height();
            BADAssSimpleDump(x, y, z, h, hconf);
        });

        BADAssBool(isBlockedPosition(x, y, z), "wrong particle position", [&] ()
        {
            const int h = solver().height(x, y);
            const double hconf = solver().confiningSurfaceEvent().height();
            BADAssSimpleDump(x, y, z, h, hconf);
        });

        for (const SOSDiffusionReaction *reaction2 : m_diffusionReactions)
        {
            if (reaction != reaction2)
            {
                bool equalX = reaction->x() == reaction2->x();
                bool equalY = reaction->y() == reaction2->y();
                bool equalZ = reaction->z() == reaction2->z();

                if (equalX && equalY && equalZ)
                {
                    BADAssBreak("two reactions are on the same spot.");
                }
            }
        }
    }
#endif

}


void LatticeDiffusion::setupInitialConditions()
{
    const double hMax = solver().confiningSurfaceEvent().height();
    const int zMin = solver().heights().min() + 1;

    const uint nLatticeParticles = solver().volume()*solver().concentration();

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
            z0 = zMin + rng.uniform()*(hMax - zMin);

        } while(solver().isBlockedPosition(x0, y0, z0) ||
                solver().isSurfaceSite(x0, y0, z0) ||
                isBlockedPosition(x0, y0, z0));

        addDiffusionReactant(x0, y0, z0);

        n++;
    }
}


void LatticeDiffusion::executeDiffusionReaction(SOSDiffusionReaction *reaction,
                                                            const uint x, const uint y, const int z)
{

    //particle has transitioned outside the regime.
    if (solver().isOutsideBox(x, y))
    {
        //derp fix concentration stuff
        removeDiffusionReactant(reaction, false);
    }


    //set x, y, z even if we remove the reaction in case the reaction is used
    //as the current reaction in the kmcsolver.
    reaction->setX(x);
    reaction->setY(y);
    reaction->setZ(z);

    int roof = 100000;
    if (z > roof)
    {
        removeDiffusionReactant(reaction, false);
        return;
    }

    if (solver().isSurfaceSite(x, y, z))
    {
        m_mutexSolver.registerHeightChange(x, y, 1);
        removeDiffusionReactant(reaction, false);

        int zAbove = z+1;
        while (isBlockedPosition(x, y, zAbove))
        {
            m_mutexSolver.registerHeightChange(x, y, 1);
            removeDiffusionReactant(x, y, zAbove, false);
//            cout << "removed " << x << " " << y << " " << zAbove << " new h = " << solver().height(x, y) << endl;
            zAbove++;
//            cout << "checking " << x << " " << y << " " << zAbove << endl;
        }

    }

}

bool LatticeDiffusion::isBlockedPosition(const uint x, const uint y, const int z) const
{
    return diffusionReaction(x, y, z) != NULL;;
}

void LatticeDiffusion::registerHeightChange(const uint x, const uint y, const int delta)
{
    //position particle on random surrounding site
    //add one to solution site heights because at this stage height(x,y) is already updated
    if (delta == -1)
    {
        const uint randomSite = rng.uniform()*solver().numberOfSurroundingSolutionSites(x, y, solver().height(x, y)+1);

        int dx, dy, dz;
        solver().getSolutionSite(x, y, solver().height(x, y)+1, dx, dy, dz, randomSite);

        const uint xNew = solver().boundaryTransform(x, dx, 0);
        const uint yNew = solver().boundaryTransform(y, dy, 1);
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

double LatticeDiffusion::depositionRate(const uint x, const uint y) const
{
    (void) x;
    (void) y;

    //Deposition is modelled as a diffusion reaction and is not explicitly treated.
    return 0;
}

