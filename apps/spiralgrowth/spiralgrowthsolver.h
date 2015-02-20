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

    void changeHeight(const uint x)
    {
        m_heights(x)++;
        m_heights.save("/tmp/heights.arma");
    }

    // KMCSolver interface
private:
    void initializeReactions();
    void updateReactions();
    void selectReaction();
    double getTotalRate();

    const uint m_length;

    vec m_heights;

    vector<TestReaction*> m_reactions;

};


class TestReaction : public Reaction
{
public:
    TestReaction(const uint i, SpiralGrowthSolver *solver) :
        Reaction(),
        m_i(i),
        m_solver(solver)
    {

    }

    virtual ~TestReaction() {}

    // Reaction interface
public:
    bool isAllowed() const
    {
        return true;
    }

    void execute()
    {
        m_solver->changeHeight(m_i);
    }

private:

    const uint m_i;

    SpiralGrowthSolver *m_solver;

};
