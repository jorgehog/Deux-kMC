#include "latticediffusion.h"

#include "../confiningsurface/confiningsurface.h"

#include "../kmcsolver/boundary/boundary.h"

#include "dissolutiondeposition.h"

#include "sosdiffusionreaction.h"

#include "concentrationboundaryreaction.h"

#include <utility>


LatticeDiffusion::LatticeDiffusion(SOSSolver &solver) :
    Diffusion(solver, "latticeDiffusion"),
    m_mutexSolver(solver)
{

}

LatticeDiffusion::~LatticeDiffusion()
{
    deleteQueuedReactions();

    for (auto &m : m_diffusionReactionsMap)
    {
        delete m.second;
    }

    //    for (SOSDiffusionReaction *reaction : m_diffusionReactionsMap)
    //    {
    //        delete reaction;
    //    }

    m_diffusionReactionsMap.clear();
}


void LatticeDiffusion::removeDiffusionReactant(SOSDiffusionReaction *reaction, bool _delete)
{
    m_diffusionReactionsMap.erase(indices(reaction));

    m_mutexSolver.removeReaction(reaction);

    //    auto &r = m_diffusionReactions;
    //    r.erase( std::remove( r.begin(), r.end(), reaction ), r.end() );

    m_mutexSolver.updateConcentrationBoundaryIfOnBoundary(reaction->x(), reaction->y());

    registerAffectedAround(reaction->x(), reaction->y(), reaction->z());

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
    const auto &res = m_diffusionReactionsMap.find(indices(x, y, z));

    if (res != m_diffusionReactionsMap.end())
    {
        return res->second;
    }

    else
    {
        return nullptr;
    }

    //    for (SOSDiffusionReaction *reaction : m_diffusionReactions)
    //    {
    //        if ((x == reaction->x()) && (y == reaction->y()) && (z == reaction->z()))
    //        {
    //            return reaction;
    //        }
    //    }

    //    return nullptr;
}

void LatticeDiffusion::clearDiffusionReactions()
{
    SOSDiffusionReaction *reaction;

    for (auto &m : m_diffusionReactionsMap)
    {
        reaction = m.second;
        m_mutexSolver.removeReaction(reaction);
        delete reaction;
    }

    m_diffusionReactionsMap.clear();
}

void LatticeDiffusion::deleteQueuedReactions()
{

    for (SOSDiffusionReaction *reaction : m_deleteQueue)
    {
        BADAssEqual(nullptr, diffusionReaction(reaction->x(), reaction->y(), reaction->z()));
        delete reaction;
    }

    m_deleteQueue.clear();
}

void LatticeDiffusion::attachToSurface(const uint x,
                                       const uint y,
                                       const int z,
                                       SOSDiffusionReaction *reaction)
{
    BADAssBool(solver().isSurfaceSite(x, y, z));

    m_mutexSolver.registerHeightChange(x, y, 1);
    removeDiffusionReactant(reaction, false);

    int zAbove = z+1;
    while (isBlockedPosition(x, y, zAbove))
    {
        m_mutexSolver.registerHeightChange(x, y, 1);
        removeDiffusionReactant(x, y, zAbove, false);
        zAbove++;
    }
}

void LatticeDiffusion::registerAffectedAround(const uint x, const uint y, const int z)
{
    const int left = solver().leftSite(x);
    const int right = solver().rightSite(x);
    const int top = solver().topSite(y);
    const int bottom = solver().bottomSite(y);

    registerAffectedAroundSingle(left, y, 0, z);
    registerAffectedAroundSingle(right, y, 0, z);
    registerAffectedAroundSingle(top, x, 1, z);
    registerAffectedAroundSingle(bottom, x, 1, z);

    SOSDiffusionReaction *r;
    if ((r = diffusionReaction(x, y, z - 1)) != nullptr)
    {
        m_mutexSolver.registerAffectedReaction(r);
    }

    if ((r = diffusionReaction(x, y, z + 1)) != nullptr)
    {
        m_mutexSolver.registerAffectedReaction(r);
    }

}

void LatticeDiffusion::registerAffectedAroundSingle(const int neighbor, const uint xi, const uint dim, const int z)
{
    if (solver().isOutsideBoxSingle(neighbor, dim))
    {
        return;
    }

    int hNeighbor;
    SOSDiffusionReaction *r;

    if (dim == 0)
    {

        //there is only a change in the rates of a surface reaction if
        //there is a change at the top of the surface (SOS)
        hNeighbor = solver().height(neighbor, xi);
        if (z == hNeighbor || z == hNeighbor + 1)
        {
            m_mutexSolver.registerAffectedReaction(&m_mutexSolver.surfaceReaction(neighbor, xi));
        }

        r = diffusionReaction(neighbor, xi, z);
        if (r != nullptr)
        {
            m_mutexSolver.registerAffectedReaction(r);
        }
    }

    else
    {
        hNeighbor = solver().height(xi, neighbor);
        if (z == hNeighbor || z == hNeighbor + 1)
        {
            m_mutexSolver.registerAffectedReaction(&m_mutexSolver.surfaceReaction(xi, neighbor));
        }

        r = diffusionReaction(xi, neighbor, z);
        if (r != nullptr)
        {
            m_mutexSolver.registerAffectedReaction(r);
        }
    }

}

void LatticeDiffusion::dumpDiffusingParticles(const uint frameNumber) const
{
    const double &h = solver().confiningSurfaceEvent().height();
    const int zMin = solver().heights().min();

    lammpswriter writer(4, "kmcdiff", "/tmp");
    writer.setSystemSize(solver().length(), solver().width(), h, 0, 0, zMin);
    writer.initializeNewFile(frameNumber);

    SOSDiffusionReaction *reaction;
    for (auto &m : m_diffusionReactionsMap)
    {
        reaction = m.second;

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

void LatticeDiffusion::moveReaction(SOSDiffusionReaction *reaction, const uint x, const uint y, const int z)
{
    m_diffusionReactionsMap.erase(indices(reaction));
    m_diffusionReactionsMap[indices(x, y, z)] = reaction;

    reaction->setX(x);
    reaction->setY(y);
    reaction->setZ(z);
}

void LatticeDiffusion::dump(const uint frameNumber) const
{
    Diffusion::dump(frameNumber);

    dumpDiffusingParticles(frameNumber);
}

SOSDiffusionReaction *LatticeDiffusion::addDiffusionReactant(const uint x, const uint y, const int z, bool setRate)
{
    //we do not add the particle if it is outside the box
    if (solver().isOutsideBox(x, y))
    {
        return nullptr;
    }

    //if we add the particle to the surface, the particle becomes part of the surface
    if (solver().isSurfaceSite(x, y, z))
    {
        m_mutexSolver.registerHeightChange(x, y, 1);
        return nullptr;
    }

    BADAssBool(!isBlockedPosition(x, y, z), "spot already taken");
    BADAssBool(!solver().isBlockedPosition(x, y, z), "spot inside surfaces");

    SOSDiffusionReaction *reaction = new SOSDiffusionReaction(m_mutexSolver, x, y, z);

    //    m_diffusionReactions.push_back(reaction);
    m_diffusionReactionsMap[indices(x, y, z)] = reaction;
//    m_diffusionReactionsMap.insert(std::make_pair(indices(x, y, z), reaction));


    m_mutexSolver.addReaction(reaction);

    m_mutexSolver.updateConcentrationBoundaryIfOnBoundary(x, y);

    reaction->setNumberOfFreePaths();

    if (setRate)
    {
        reaction->calculateRate();
    }

    registerAffectedAround(x, y, z);

    return reaction;
}

void LatticeDiffusion::execute()
{
    deleteQueuedReactions();

#ifndef NDEBUG
    SOSDiffusionReaction *reaction;
    for (auto &m : m_diffusionReactionsMap)
    {
        reaction = m.second;

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

        SOSDiffusionReaction *reaction2;
        for (auto &m2 : m_diffusionReactionsMap)
        {
            reaction2 = m2.second;
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

void LatticeDiffusion::reset()
{
    //    for (SOSDiffusionReaction *r : m_diffusionReactions)
    //    {
    //        BADAssEqual(m_diffusionReactionsMap[r->x()][r->y()][r->z()], r);
    //    }
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

        addDiffusionReactant(x0, y0, z0, false);

        n++;
    }
}


void LatticeDiffusion::executeDiffusionReaction(SOSDiffusionReaction *reaction,
                                                const int x, const int y, const int z)
{
    const uint xOld = reaction->x();
    const uint yOld = reaction->y();
    const int zOld = reaction->z();

    registerAffectedAround(xOld, yOld, zOld);

    //particle has transitioned outside the regime.
    if (solver().isOutsideBox(x, y))
    {
        removeDiffusionReactant(reaction, false);
        return;
    }

    const uint ux = (uint)x;
    const uint uy = (uint)y;

    registerAffectedAround(ux, uy, z);

    moveReaction(reaction, ux, uy, z);

    if (solver().isSurfaceSite(ux, uy, z))
    {
        attachToSurface(ux, uy, z, reaction);
    }

    if (!(xOld == ux && yOld == uy))
    {
        m_mutexSolver.updateConcentrationBoundaryIfOnBoundary(xOld, yOld);
        m_mutexSolver.updateConcentrationBoundaryIfOnBoundary(ux, uy);
    }

    m_mutexSolver.registerAffectedReaction(reaction);

}

void LatticeDiffusion::executeConcentrationBoundaryReaction(ConcentrationBoundaryReaction *reaction)
{
    const uint n =reaction->freeBoundarySites();
    const uint nChosen = rng.uniform()*n;

    uint xi;
    int z;
    reaction->getFreeBoundarSite(nChosen, xi, z);

    if (reaction->dim() == 0)
    {
        addDiffusionReactant(reaction->location(), xi, z);
    }

    else
    {
        addDiffusionReactant(xi, reaction->location(), z);
    }
}

bool LatticeDiffusion::isBlockedPosition(const uint x, const uint y, const int z) const
{
    return diffusionReaction(x, y, z) != nullptr;
}

void LatticeDiffusion::registerHeightChange(const uint x, const uint y, const int delta)
{
    //position particle on random surrounding site
    //add one to solution site heights because at this stage height(x,y) is already updated
    if (delta == -1)
    {
        const uint nSites = solver().numberOfSurroundingSolutionSites(x, y, solver().height(x, y)+1);
        const uint randomSite = rng.uniform()*nSites;

        BADAss(nSites, !=, 0u);

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

        registerAffectedAround(x, y, solver().height(x, y) + 1);
    }

    //this can occur if the surface is changed so that it connects to a particle in the solution.
    else
    {
        const int zSurface = solver().height(x, y) + 1;
        SOSDiffusionReaction *r = diffusionReaction(x, y, zSurface);

        registerAffectedAround(x, y, solver().height(x, y));

        if (r != nullptr)
        {
            attachToSurface(x, y, zSurface, r);
        }
    }

}

double LatticeDiffusion::depositionRate(const uint x, const uint y) const
{
    (void) x;
    (void) y;

    //Deposition is modelled as a diffusion reaction and is not explicitly treated.
    return 0;
}

uint LatticeDiffusion::dissolutionPaths(const uint x, const uint y) const
{
    return solver().numberOfSurroundingSolutionSites(x, y);
}


indices::indices(const SOSDiffusionReaction *r) :
    m_x(r->x()),
    m_y(r->y()),
    m_z(r->z())
{

}
