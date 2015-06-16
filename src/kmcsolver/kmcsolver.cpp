#include "kmcsolver.h"

#include "reaction.h"

#include "boundary/boundary.h"


using namespace kMC;


KMCSolver::KMCSolver(vector<const Boundary *> boundaries) :
    ignis::LatticeEvent("KMCSolver")
{
    m_boundaries = boundaries;
}


KMCSolver::~KMCSolver()
{
    m_boundaries.clear();
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

        m_cumsumRates.at(i) = m_totalRate;
    }
}

double KMCSolver::getRandomLogNumber() const
{
    return -std::log(rng.uniform());
}

void KMCSolver::getCumsumAndTotalRate()
{
    m_totalRate = 0;

    double rate;
    Reaction *reaction;

    for (uint i = 0; i < numberOfReactions(); ++i)
    {
        reaction = getReaction(i);

        if (reaction->isAllowed())
        {
            rate = reaction->rate();
            BADAssClose(reaction->rate(), reaction->rateExpression(), 1E-5);

            m_totalRate += rate;
        }
        else
        {
            rate = 0;
        }

        m_cumsumRates.at(i) = m_totalRate;
    }
}


void KMCSolver::updateTime()
{
    setCurrentTimeStep(m_nextRandomLogNumber/m_totalRate);

    BADAss(currentTimeStep(), >, 0, "timestep should be posiive.");

    m_currentTime += currentTimeStep();
}


void KMCSolver::initialize()
{
    m_currentTime = 0;
    m_nextRandomLogNumber = getRandomLogNumber();

    initializeSolver();

    initializeReactions();

    updateTime();
}


void KMCSolver::execute()
{
    double R = rng.uniform()*m_totalRate;

    uint choice = binarySearchForInterval(R, m_cumsumRates);

    m_selectedReaction = getReaction(choice);

    m_nextRandomLogNumber = getRandomLogNumber();

    BADAssBool(selectedReaction()->isAllowed());
}


void KMCSolver::reset()
{
    m_selectedReaction->executeAndUpdate();

    getCumsumAndTotalRate();

    updateTime();
}

