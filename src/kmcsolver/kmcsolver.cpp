#include "kmcsolver.h"

#include "reaction.h"

#include "boundary/boundary.h"


using namespace kMC;


KMCSolver::KMCSolver() :
    ignis::LatticeEvent("KMCSolver")
{

}


KMCSolver::~KMCSolver()
{
    //KMCSolver is not responsible for deletion
    m_reactions.clear();
    m_affectedReactions.clear();
    m_boundaries.clear();
}


void KMCSolver::setCurrentTimeStep(double currentTimeStep)
{
    m_currentTimeStep = currentTimeStep;
}

void KMCSolver::resizeCumsumRates()
{
    if (m_cumsumRates.size() < numberOfReactions())
    {
        m_cumsumRates.resize(numberOfReactions());
    }
}

void KMCSolver::initializeReactions()
{
    updateAffectedReactions();

//    m_totalRate = 0;

//    m_cumsumRates.resize(numberOfReactions());

//    double rate;
    Reaction *reaction;

    for (uint i = 0; i < numberOfReactions(); ++i)
    {
        reaction = getReaction(i);

        if (reaction->isAllowed())
        {
            reaction->calculateRate();
//            rate = reaction->rate();

//            m_totalRate += rate;
        }
//        else
//        {
//            rate = 0;
//        }

//        m_cumsumRates.at(i) = m_totalRate;
    }
}

void KMCSolver::updateAffectedReactions()
{
    for (Reaction *r : m_affectedReactions)
    {
        r->affectedUpdateRule();
    }

    m_affectedReactions.clear();
}

void KMCSolver::addReaction(Reaction *reaction)
{
    m_reactions.push_back(reaction);
}

void KMCSolver::removeReaction(const uint n)
{
    auto it = m_reactions.begin() + n;
    Reaction *reaction = *it;

    m_reactions.erase(it);

    removeIfAffected(reaction);
}

void KMCSolver::removeReaction(Reaction *reaction)
{
    auto &r = m_reactions;
    r.erase( std::remove( r.begin(), r.end(), reaction ), r.end() );

    removeIfAffected(reaction);
}

void KMCSolver::removeIfAffected(Reaction *reaction)
{
    auto it = std::find(m_affectedReactions.begin(), m_affectedReactions.end(), reaction);

    if (it != m_affectedReactions.end())
    {
        m_affectedReactions.erase(it);
    }
}

void KMCSolver::recalculateAllRates()
{
    for (Reaction *reaction : m_reactions)
    {
        if (reaction->isAllowed())
        {
            reaction->calculateRate();
        }
    }
}

double KMCSolver::timeUnit() const
{
    return 1.0;
}

double KMCSolver::getRandomLogNumber() const
{
    return -std::log(rng.uniform());
}

void KMCSolver::getCumsumAndTotalRate()
{
    updateAffectedReactions();

    m_totalRate = 0;

    resizeCumsumRates();

    uint i = 0;
    for (const Reaction *reaction : m_reactions)
    {
        if (reaction->isAllowed())
        {
            m_totalRate += reaction->rate();
        }

        m_cumsumRates[i] = m_totalRate;

        i++;
    }

}


void KMCSolver::updateTime()
{
    setCurrentTimeStep(m_nextRandomLogNumber/(m_totalRate*timeUnit()));

    BADAss(currentTimeStep(), >, 0, "timestep should be posiive.");

    m_currentTime += currentTimeStep();
}


void KMCSolver::initialize()
{
    BADAssBool(!m_boundaries.empty(), "boundaries are not set.");

    m_currentTime = 0;
    m_nextRandomLogNumber = getRandomLogNumber();

    initializeReactions();

//    updateTime();
}


void KMCSolver::execute()
{
    m_nextRandomLogNumber = getRandomLogNumber();

    getCumsumAndTotalRate();

    updateTime();

    BADAss(m_cumsumRates.size(), >=, numberOfReactions());

    uint choice = chooseFromTotalRate(m_cumsumRates, m_totalRate, numberOfReactions());

    m_selectedReaction = getReaction(choice);

    BADAssBool(selectedReaction()->isAllowed());
}


void KMCSolver::reset()
{
    m_selectedReaction->executeAndUpdate();
}

