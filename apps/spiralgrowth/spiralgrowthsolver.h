#pragma once

#include <kMC>
#include <armadillo>

using namespace kMC;
using namespace arma;

class TestReaction;

class SpiralGrowthSolver : public KMCSolver<uint>
{
public:
    SpiralGrowthSolver(const uint length);
    ~SpiralGrowthSolver();

    // KMCSolver interface
private:
    void initializeReactions();
    void updateReactions();
    void selectReaction();
    void executeReaction(uint &reaction);
    double getTotalRate();

    const uint m_length;

    vec m_heights;

    uint m_reaction;

};
