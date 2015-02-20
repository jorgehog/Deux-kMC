#pragma once

#include <ignis/include/ignis.h>

#include "kmcassets.h"

#include "RNG/rng.h"

namespace kMC
{

class Reaction;

class KMCSolver : public LatticeEvent
{
public:
    KMCSolver();
    ~KMCSolver();

    double currentTimeStep() const;

private:

    Reaction *m_selectedReaction;

    double m_currentTimeStep;
    double m_currentTime;

    virtual void initializeReactions() = 0;

    virtual void updateReactions() = 0;

    virtual void selectReaction() = 0;

    virtual double getTotalRate() = 0;


    void updateTime();

    void setCurrentTimeStep(double currentTimeStep);

protected:

    inline void setSelectedReaction(Reaction *reaction);

    // Event interface
public:
    virtual void initialize() override final;

    virtual void execute() override final;

    virtual void reset() override final;

};

void KMCSolver::setSelectedReaction(Reaction *reaction)
{
    m_selectedReaction = reaction;
}

} //End of namespace kMC
