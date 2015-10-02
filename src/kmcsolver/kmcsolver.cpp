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

void KMCSolver::initializeReactions()
{
    updateAffectedReactions();

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
    m_reactions.erase(m_reactions.begin() + n);
}

void KMCSolver::removeReaction(Reaction *reaction)
{
    auto &r = m_reactions;
    r.erase( std::remove( r.begin(), r.end(), reaction ), r.end() );
}

double KMCSolver::getRandomLogNumber() const
{
    return -std::log(rng.uniform());
}

void KMCSolver::getCumsumAndTotalRate()
{
    updateAffectedReactions();

    m_totalRate = 0;

    double rate;
    Reaction *reaction;

    if (m_cumsumRates.size() != numberOfReactions())
    {
        m_cumsumRates.resize(numberOfReactions());
    }

    for (uint i = 0; i < numberOfReactions(); ++i)
    {
        reaction = getReaction(i);

        if (reaction->isAllowed())
        {
            rate = reaction->rate();
            BADAssClose(reaction->rate(),
                        reaction->rateExpression(),
                        1E-5,
                        "Reaction rate inconsistency.",
                        [&] ()
            {
                double r1 = reaction->rate();
                double r2 = reaction->rateExpression();

                BADAssSimpleDump(r1, r2);
            });

            m_totalRate += rate;
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
    BADAssBool(!m_boundaries.empty(), "boundaries are not set.");

    m_currentTime = 0;
    m_nextRandomLogNumber = getRandomLogNumber();

    initializeReactions();

    updateTime();
}


void KMCSolver::execute()
{
    uint choice = chooseFromTotalRate(m_cumsumRates, m_totalRate);

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

