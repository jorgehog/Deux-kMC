#pragma once

#include <ignis/include/ignis.h>

#include "kmcassets.h"

#include "RNG/rng.h"

namespace kMC
{

class Reaction;

class KMCSolver : public ignis::LatticeEvent
{
public:
    KMCSolver();
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

    virtual uint numberOfReactions() const = 0;

    virtual Reaction *getReaction(const uint n) const = 0;

private:

    Reaction *m_selectedReaction;

    double m_currentTimeStep;
    double m_currentTime;

    double m_totalRate;

    vector<double> m_cumsumRates;

    void updateTime();

    void setCurrentTimeStep(double currentTimeStep);

    void initializeReactions();

    void getCumsumAndTotalRate();

    // Event interface
public:
    virtual void initialize() override final;

    virtual void execute() override final;

    virtual void reset() override final;

};

} //End of namespace kMC
