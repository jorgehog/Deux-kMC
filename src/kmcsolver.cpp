#include "kmcsolver.h"

#include "reaction.h"


using namespace kMC;

KMCSolver::KMCSolver()
{

}

KMCSolver::~KMCSolver()
{

}

double KMCSolver::currentTimeStep() const
{
    return m_currentTimeStep;
}

void KMCSolver::setCurrentTimeStep(double currentTimeStep)
{
    m_currentTimeStep = currentTimeStep;
}

void KMCSolver::updateTime()
{
    setCurrentTimeStep(-std::log(rng.uniform())/getTotalRate());

    BADAss(currentTimeStep(), >, 0, "timestep should be posiive.");

    m_currentTime += currentTimeStep();
}

void KMCSolver::initialize()
{
    initializeReactions();

    m_currentTime = 0;

    updateTime();
}

void KMCSolver::execute()
{
    selectReaction();
}

void KMCSolver::reset()
{
    m_selectedReaction->execute();

    updateReactions();

    updateTime();
}

