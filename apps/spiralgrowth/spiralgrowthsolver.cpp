#include "spiralgrowthsolver.h"

SpiralGrowthSolver::SpiralGrowthSolver(const uint length) :
    KMCSolver(),
    m_length(length),
    m_heights(length),
    m_reactions(length)
{

}

SpiralGrowthSolver::~SpiralGrowthSolver()
{
    for (uint x = 0; x < m_length; ++x)
    {
        delete m_reactions.at(x);
    }

    m_reactions.clear();
}



void SpiralGrowthSolver::initializeReactions()
{
    for (uint x = 0; x < m_length; ++x)
    {
        m_reactions.at(x) = new TestReaction(x, this);
    }
}

void SpiralGrowthSolver::updateReactions()
{

}

void SpiralGrowthSolver::selectReaction()
{
    vec ratesumvec = cumsum(cumsum(ones(m_length)));

    double R = rng.uniform()*ratesumvec(m_length-1);

    uint choice = binarySearchForInterval(R, ratesumvec);

    setSelectedReaction(m_reactions.at(choice));
}

double SpiralGrowthSolver::getTotalRate()
{
    return m_heights.size();
}

