#include "particlenumberconservator.h"

#include "sossolver.h"
#include "Events/diffusion/diffusion.h"

using namespace kMC;

ParticleNumberConservator::ParticleNumberConservator(SOSSolver &solver) :
    m_solver(solver)
{

}



void ParticleNumberConservator::initializeObserver(const Subjects &subject)
{
    (void) subject;

    BADAssBool(m_solver.diffusionEvent().hasDiscreteParticles());

    m_targetN = m_solver.diffusionEvent().numberOfParticles();

}

void ParticleNumberConservator::notifyObserver(const Subjects &subject)
{
    (void) subject;

    int particleDifference = (int)m_solver.diffusionEvent().numberOfParticles() - (int)m_targetN;

    if (particleDifference < 0)
    {
        for (int n = particleDifference; n < 0; ++n)
        {
            m_solver.diffusionEvent().insertRandomParticle();
        }
    }

    else
    {
        for (int n = 0; n < particleDifference; ++n)
        {
            m_solver.diffusionEvent().removeRandomParticle();
        }
    }

}
