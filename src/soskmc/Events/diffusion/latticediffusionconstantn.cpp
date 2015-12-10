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

void LatticeDiffusionConstantN::notifyObserver(const Subjects &subject)
{
    LatticeDiffusion::notifyObserver(subject);

    int particleDifference = numberOfDiffusionReactions() - m_targetNParticles;

    if (particleDifference < 0)
    {
        addRandomParticles(-particleDifference, true);
    }

    else
    {
        removeRandomParticles(particleDifference);
    }
}

void LatticeDiffusionConstantN::initializeObserver(const Subjects &subject)
{
    LatticeDiffusion::initializeObserver(subject);

    m_targetNParticles = numberOfDiffusionReactions();

}
