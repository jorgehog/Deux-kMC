#include "latticediffusionconstantn.h"

#include "../../sossolver.h"
#include "../confiningsurface/confiningsurface.h"
#include "../../sosdiffusionreaction.h"
#include "../../dissolutiondeposition.h"

LatticeDiffusionConstantN::LatticeDiffusionConstantN(SOSSolver &solver) :
    Diffusion(solver, "LatticeDiffusionConstantN", "", true, true),
    LatticeDiffusion(solver)
{

}



void LatticeDiffusionConstantN::notifyObserver()
{
    LatticeDiffusion::notifyObserver();

    int surfaceContactHeight = solver().confiningSurfaceEvent().height() - 1;

    vector<SOSDiffusionReaction *> blockedReactions;

    SOSDiffusionReaction *r;
    for (auto &m : m_diffusionReactionsMap)
    {
        r = m.second;

        if (solver().isBlockedPosition(r->x(), r->y(), r->z()))
        {
            blockedReactions.push_back(r);
        }

        if (r->z() == surfaceContactHeight)
        {
            solver().registerAffectedReaction(r);
        }
    }

    for (SOSDiffusionReaction *r : blockedReactions)
    {
        removeDiffusionReactant(r);
    }

    int particleDifference = numberOfDiffusionReactions() - m_targetNParticles;

    if (particleDifference < 0)
    {
        addRandomParticles(-particleDifference, true);
    }

    else
    {
        removeRandomParticles(particleDifference);
    }

    //very slow
    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            if (solver().height(x, y) == surfaceContactHeight)
            {
                solver().registerAffectedReaction(&solver().surfaceReaction(x, y));
            }
        }
    }



}

void LatticeDiffusionConstantN::initializeObserver()
{
    LatticeDiffusion::initializeObserver();

    m_targetNParticles = numberOfDiffusionReactions();
    m_targetSeparation = solver().confiningSurfaceEvent().height() - solver().averageHeight();

}
