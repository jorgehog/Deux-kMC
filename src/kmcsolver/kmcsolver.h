#pragma once

#include <ignis/include/ignis.h>

#include "kmcassets.h"

#include "RNG/rng.h"

namespace kMC
{

class Reaction;
class Boundary;

class KMCSolver : public ignis::LatticeEvent
{
public:
    KMCSolver(vector<const Boundary*> boundaries);
    virtual ~KMCSolver();

    const double &currentTimeStep() const
    {
        return m_currentTimeStep;
    }

    const double &currentTime() const
    {
        return m_currentTime;
    }

    const Reaction *selectedReaction() const
    {
        return m_selectedReaction;
    }

    const Boundary *boundary(const uint i) const
    {
        return m_boundaries.at(i);
    }

    void initializeReactions();

    Reaction *getReaction(const uint n) const
    {
        return m_reactions.at(n);
    }

    const double &nextRandomLogNumber() const
    {
        return m_nextRandomLogNumber;
    }

    void addReaction(Reaction *reaction);

    void removeReaction(const uint n);

    void removeReaction(Reaction *reaction);

    uint numberOfReactions() const
    {
        return m_reactions.size();
    }

private:
    vector<Reaction*> m_reactions;

    vector<const Boundary*> m_boundaries;

    Reaction *m_selectedReaction;

    double m_currentTimeStep;
    double m_currentTime;

    double m_totalRate;

    double m_nextRandomLogNumber;

    vector<double> m_cumsumRates;

    double getRandomLogNumber() const;

    void updateTime();

    void setCurrentTimeStep(double currentTimeStep);

    void getCumsumAndTotalRate();

private:
    virtual void initializeSolver() = 0;

    // Event interface
public:
    virtual void initialize() override final;

    virtual void execute() override final;

    virtual void reset() override final;

};

} //End of namespace kMC
