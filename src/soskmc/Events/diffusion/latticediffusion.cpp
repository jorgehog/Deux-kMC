#include "latticediffusion.h"

#include "../confiningsurface/confiningsurface.h"

#include "../kmcsolver/boundary/boundary.h"

#include "dissolutiondeposition.h"

#include "sosdiffusionreaction.h"

#include "concentrationboundaryreaction.h"

#include <utility>


LatticeDiffusion::LatticeDiffusion(SOSSolver &solver) :
    Diffusion(solver, "latticeDiffusion", "", true, true)
{

}

LatticeDiffusion::~LatticeDiffusion()
{
    deleteQueuedReactions();

    for (auto &m : m_diffusionReactionsMap)
    {
        delete m.second;
    }

    m_diffusionReactionsMap.clear();
}


void LatticeDiffusion::removeDiffusionReactant(SOSDiffusionReaction *reaction, bool _delete)
{
    BADAss(std::find(m_deleteQueue.begin(), m_deleteQueue.end(), reaction), ==, m_deleteQueue.end());
    BADAssBool(isBlockedPosition(reaction->x(), reaction->y(), reaction->z()), "supposed to remove non existing particle", [&reaction] ()
    {
        BADAssSimpleDump(reaction->x(), reaction->y(), reaction->z());
    });

    m_diffusionReactionsMap.erase(indices(reaction));

    solver().removeReaction(reaction);

    solver().updateConcentrationBoundaryIfOnBoundary(reaction->x(), reaction->y());

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

void LatticeDiffusion::removeDiffusionReactant(const int x, const int y, const int z, bool _delete)
{
    removeDiffusionReactant(diffusionReaction(x, y, z), _delete);
}

void LatticeDiffusion::removeDiffusionReactant(const uint index, bool _delete)
{
    auto it = m_diffusionReactionsMap.begin();

    std::advance(it, index);

    removeDiffusionReactant((*it).second, _delete);
}

void LatticeDiffusion::removeRandomParticles(const uint N, bool _delete)
{
    uint n = 0;
    while (n < N)
    {
        const uint which = rng.uniform()*numberOfDiffusionReactions();

        removeDiffusionReactant(which, _delete);

        n++;
    }
}

void LatticeDiffusion::addRandomParticles(const uint N, bool setRate)
{
    const double hMax = solver().confiningSurfaceEvent().height();
    const int zMin = solver().heights().min() + 1;

    uint x0;
    uint y0;
    int z0;

    uint n = 0;
    while (n < N)
    {
        do
        {
            x0 = rng.uniform()*solver().length();
            y0 = rng.uniform()*solver().width();
            z0 = floor(zMin + rng.uniform()*(hMax - zMin));

        } while(solver().isBlockedPosition(x0, y0, z0) ||
                isBlockedPosition(x0, y0, z0) ||
                solver().isSurfaceSite(x0, y0, z0));

        addDiffusionReactant(x0, y0, z0, setRate);

        BADAssBool(!solver().isSurfaceSite(x0, y0, z0));

        n++;
    }
}

SOSDiffusionReaction *LatticeDiffusion::diffusionReaction(const int x, const int y, const int z) const
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
}

void LatticeDiffusion::clearDiffusionReactions()
{
    SOSDiffusionReaction *reaction;

    for (auto &m : m_diffusionReactionsMap)
    {
        reaction = m.second;
        solver().removeReaction(reaction);
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
                                       const int z)
{
    int zAbove = z;
    while (isBlockedPosition(x, y, zAbove))
    {
        BADAssBool(solver().isSurfaceSite(x, y, zAbove));
        BADAssBool(isBlockedPosition(x, y, zAbove));

        removeDiffusionReactant(x, y, zAbove, false);
        solver().registerHeightChange(x, y, 1);
        zAbove++;
    }
}

void LatticeDiffusion::registerAffectedAround(const uint x, const uint y, const int z)
{
    const int left = solver().leftSite(x, y, z);
    const int right = solver().rightSite(x, y, z);
    const int top = solver().topSite(x, y, z);
    const int bottom = solver().bottomSite(x, y, z);

    registerAffectedAroundSingle(left, y, 0, z);
    registerAffectedAroundSingle(right, y, 0, z);
    registerAffectedAroundSingle(top, x, 1, z);
    registerAffectedAroundSingle(bottom, x, 1, z);

    SOSDiffusionReaction *r;
    if ((r = diffusionReaction(x, y, z - 1)) != nullptr)
    {
        solver().registerAffectedReaction(r);
    }

    if ((r = diffusionReaction(x, y, z + 1)) != nullptr)
    {
        solver().registerAffectedReaction(r);
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
            solver().registerAffectedReaction(&solver().surfaceReaction(neighbor, xi));
        }

        r = diffusionReaction(neighbor, xi, z);
        if (r != nullptr)
        {
            solver().registerAffectedReaction(r);
        }
    }

    else
    {
        hNeighbor = solver().height(xi, neighbor);
        if (z == hNeighbor || z == hNeighbor + 1)
        {
            solver().registerAffectedReaction(&solver().surfaceReaction(xi, neighbor));
        }

        r = diffusionReaction(xi, neighbor, z);
        if (r != nullptr)
        {
            solver().registerAffectedReaction(r);
        }
    }

}

void LatticeDiffusion::dumpDiffusingParticles(const uint frameNumber, const string path) const
{
    const double &h = solver().confiningSurfaceEvent().height();
    const int zMin = solver().heights().min();

    lammpswriter writer(4, "kmcdiff", path);
    writer.setSystemSize(solver().length(), solver().width(), h, 0, 0, zMin);
    writer.initializeNewFile(frameNumber);

    if (m_diffusionReactionsMap.empty())
    {
        writer << 3 << 0 << 0 << solver().height(0, 0);
    }

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

vector<SOSDiffusionReaction *> LatticeDiffusion::particlesSurrounding(const uint x, const uint y, const int z) const
{
    int left = solver().leftSite(x, y, z);
    int right = solver().rightSite(x, y, z);
    int top = solver().topSite(x, y, z);
    int bottom = solver().bottomSite(x, y, z);

    int below = z - 1;
    int above = z + 1;

    vector<SOSDiffusionReaction *> particles;

    SOSDiffusionReaction *r;

    if ((r = diffusionReaction(left, y, z)) != nullptr)
    {
        particles.push_back(r);
    }

    if ((r = diffusionReaction(right, y, z)) != nullptr)
    {
        particles.push_back(r);
    }

    if ((r = diffusionReaction(x, top, z)) != nullptr)
    {
        particles.push_back(r);
    }

    if ((r = diffusionReaction(x, bottom, z)) != nullptr)
    {
        particles.push_back(r);
    }

    if ((r = diffusionReaction(x, y, below)) != nullptr)
    {
        particles.push_back(r);
    }

    if ((r = diffusionReaction(x, y, above)) != nullptr)
    {
        particles.push_back(r);
    }

    return particles;


}

void LatticeDiffusion::dump(const uint frameNumber, const string path) const
{
    Diffusion::dump(frameNumber, path);

    dumpDiffusingParticles(frameNumber, path);
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
        solver().registerHeightChange(x, y, 1);
        return nullptr;
    }

    BADAssBool(!isBlockedPosition(x, y, z), "spot already taken");
    BADAssBool(!solver().isBlockedPosition(x, y, z), "spot inside surfaces");

    SOSDiffusionReaction *reaction = new SOSDiffusionReaction(solver(), x, y, z);

    m_diffusionReactionsMap[indices(x, y, z)] = reaction;

    solver().addReaction(reaction);

    solver().updateConcentrationBoundaryIfOnBoundary(x, y);

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
    setValue(numberOfDiffusionReactions());

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

void LatticeDiffusion::setupInitialConditions()
{

    //already initialized particles
    if (numberOfDiffusionReactions() != 0)
    {
        return;
    }

    //subtract area from volume since we do not initiate surface particles
    //add a random number such that if we get 3.3 particles, there is a 0.3 chance to get 3 + 1.
    const uint nLatticeParticles = solver().freeVolume()*solver().concentration() + rng.uniform();

    addRandomParticles(nLatticeParticles);
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
        attachToSurface(ux, uy, z);
    }

    else
    {
        solver().registerAffectedReaction(reaction);
    }


    if (!(xOld == ux && yOld == uy))
    {
        solver().updateConcentrationBoundaryIfOnBoundary(xOld, yOld);
        solver().updateConcentrationBoundaryIfOnBoundary(ux, uy);
    }
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

double LatticeDiffusion::concentration() const
{
    return m_diffusionReactionsMap.size()/solver().freeVolume();
}

void LatticeDiffusion::registerHeightChange(const uint x,
                                            const uint y,
                                            const int value,
                                            std::vector<DissolutionDeposition *> &affectedSurfaceReactions,
                                            const uint n)
{
    (void) affectedSurfaceReactions;
    (void) n;

    //position particle on random surrounding site
    //add one to solution site heights because at this stage height(x,y) is already updated
    if (value == -1)
    {
        const int h = solver().height(x, y)+1;

        int dx, dy, dz;
        solver().getRandomSolutionSite(x, y, h, dx, dy, dz);

        const uint xNew = solver().boundaryTransform(x, y, h, dx, 0);
        const uint yNew = solver().boundaryTransform(x, y, h, dy, 1);
        const int zNew = h + dz;

        BADAssBool(!solver().isSurfaceSite(xNew, yNew, zNew), "adding diff reaction to surface.", [&] ()
        {
            int hp = solver().height(x, y);
            int hNew = solver().height(xNew, yNew);
            BADAssSimpleDump(x, y, hp, xNew, yNew, zNew, hNew);
        });

        addDiffusionReactant(xNew, yNew, zNew);

        registerAffectedAround(x, y, h);
    }

    //this can occur if the surface is changed so that it connects to a particle in the solution.
    else
    {
//        const int zSurface = solver().height(x, y) + 1;
//        SOSDiffusionReaction *r = diffusionReaction(x, y, zSurface);

        registerAffectedAround(x, y, solver().height(x, y));

//        if (r != nullptr)
//        {
//            attachToSurface(x, y, zSurface, r);
//        }
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
