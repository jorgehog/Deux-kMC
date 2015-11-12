#include "latticediffusionrescaling.h"

#include "../../sossolver.h"
#include "../confiningsurface/confiningsurface.h"
#include "../../sosdiffusionreaction.h"
#include "../../dissolutiondeposition.h"

LatticeDiffusionRescaling::LatticeDiffusionRescaling(SOSSolver &solver) :
    Diffusion(solver, "LatticeDiffusionRescaling", "", true, true),
    LatticeDiffusion(solver)
{

}



void LatticeDiffusionRescaling::registerHeightChange(const uint x, const uint y, const int value, std::vector<DissolutionDeposition *> &affectedSurfaceReactions, const uint nAffectedSurfaceReactions)
{
    LatticeDiffusion::registerHeightChange(x, y, value, affectedSurfaceReactions, nAffectedSurfaceReactions);

    double newHeight = solver().averageHeight() + m_targetSeparation;
    int surfaceContactHeight = (int)newHeight - 1;

    //set height to wanted separation and remove blocked particles
    solver().confiningSurfaceEvent().setHeight(newHeight);

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

void LatticeDiffusionRescaling::setupInitialConditions()
{
    LatticeDiffusion::setupInitialConditions();

    m_targetNParticles = numberOfDiffusionReactions();
    m_targetSeparation = solver().confiningSurfaceEvent().height() - solver().averageHeight();

}
