#include "spiralgrowthsolver.h"

SpiralGrowthSolver::SpiralGrowthSolver(const uint length) :
    KMCSolver<uint>(),
    m_length(length),
    m_heights(length)
{

}

SpiralGrowthSolver::~SpiralGrowthSolver()
{

}



void SpiralGrowthSolver::initializeReactions()
{

}

void SpiralGrowthSolver::updateReactions()
{

}

void SpiralGrowthSolver::selectReaction()
{
    vec ratesumvec = cumsum(cumsum(ones(m_length)));

    ratesumvec = ratesumvec%ratesumvec;

    double R = rng.uniform()*ratesumvec(m_length-1);

    m_reaction = binarySearchForInterval(R, ratesumvec);

    setSelectedReaction(m_reaction);
}

void SpiralGrowthSolver::executeReaction(uint &reaction)
{
    m_heights(reaction)++;
    m_heights.save("/tmp/heights.arma");
}

double SpiralGrowthSolver::getTotalRate()
{
    return m_heights.size();
}

