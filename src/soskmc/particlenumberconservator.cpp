#include "particlenumberconservator.h"

#include "sossolver.h"
#include "Events/diffusion/diffusion.h"

using namespace kMC;

ParticleNumberConservator::ParticleNumberConservator(SOSSolver &solver) :
    LatticeEvent("NConservator", "", true),
    m_solver(solver)
{
    BADAssBool(solver.diffusionEventIsSet(), "Diffusion event need to be set before this object can be constructed.");

    setDependency(solver.diffusionEvent());
}

void ParticleNumberConservator::initialize()
{
    BADAssBool(m_solver.diffusionEvent().hasDiscreteParticles());

    m_targetN = m_solver.diffusionEvent().numberOfParticles();
}

void ParticleNumberConservator::execute()
{
    setValue(solver().diffusionEvent().numberOfParticles());
}

void ParticleNumberConservator::reset()
{
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
