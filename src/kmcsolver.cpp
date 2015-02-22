#include "kmcsolver.h"

#include "reaction.h"

using namespace kMC;


KMCSolver::KMCSolver() :
    ignis::Event<uint>("KMCSolver")
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

void KMCSolver::initializeReactions()
{
    m_totalRate = 0;
    m_cumsumRates.resize(numberOfReactions());

    double rate;
    Reaction *reaction;

    for (uint i = 0; i < numberOfReactions(); ++i)
    {
        reaction = getReaction(i);

        if (reaction->isAllowed())
        {
            reaction->calculateRate();
            rate = reaction->rate();

            m_totalRate += rate;
        }
        else
        {
            rate = 0;
        }

        m_cumsumRates[i] = m_totalRate;
    }
}


void KMCSolver::updateTime()
{
    setCurrentTimeStep(-std::log(rng.uniform())/m_totalRate);

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
    double R = rng.uniform()*m_totalRate;

    uint choice = binarySearchForInterval(R, m_cumsumRates);

    m_selectedReaction = getReaction(choice);

}


void KMCSolver::reset()
{
    m_selectedReaction->executeAndUpdate();

//    getCumsumAndTotalRate();

    updateTime();
}

