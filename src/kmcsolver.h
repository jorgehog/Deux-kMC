#pragma once

#include <ignis/include/ignis.h>

#include "kmcassets.h"

#include "RNG/rng.h"

#include "reaction.h"

namespace kMC
{

template<typename pT>
class KMCSolver : public ignis::Event<pT>
{
public:
    KMCSolver();
    virtual ~KMCSolver();

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

template<typename pT>
void KMCSolver<pT>::setSelectedReaction(Reaction *reaction)
{
    m_selectedReaction = reaction;
}


template<typename pT>
KMCSolver<pT>::KMCSolver()
{

}

template<typename pT>
KMCSolver<pT>::~KMCSolver()
{

}

template<typename pT>
double KMCSolver<pT>::currentTimeStep() const
{
    return m_currentTimeStep;
}

template<typename pT>
void KMCSolver<pT>::setCurrentTimeStep(double currentTimeStep)
{
    m_currentTimeStep = currentTimeStep;
}

template<typename pT>
void KMCSolver<pT>::updateTime()
{
    setCurrentTimeStep(-std::log(rng.uniform())/getTotalRate());

    BADAss(currentTimeStep(), >, 0, "timestep should be posiive.");

    m_currentTime += currentTimeStep();
}

template<typename pT>
void KMCSolver<pT>::initialize()
{
    initializeReactions();

    m_currentTime = 0;

    updateTime();
}

template<typename pT>
void KMCSolver<pT>::execute()
{
    selectReaction();
}

template<typename pT>
void KMCSolver<pT>::reset()
{
    m_selectedReaction->execute();

    updateReactions();

    updateTime();
}

} //End of namespace kMC
