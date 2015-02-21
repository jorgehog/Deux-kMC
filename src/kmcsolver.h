#pragma once

#include <ignis/include/ignis.h>

#include "kmcassets.h"

#include "RNG/rng.h"

namespace kMC
{

template<class rT, typename pT=uint>
class KMCSolver : public ignis::Event<pT>
{
public:
    KMCSolver();
    virtual ~KMCSolver();

    double currentTimeStep() const;

private:

    rT *m_selectedReaction;

    double m_currentTimeStep;
    double m_currentTime;

    virtual void initializeReactions() = 0;

    virtual void updateReactions() = 0;

    virtual void selectReaction() = 0;

    virtual void executeReaction(rT &reaction) = 0;

    virtual double getTotalRate() = 0;


    void updateTime();

    void setCurrentTimeStep(double currentTimeStep);

protected:

    inline void setSelectedReaction(rT &reaction);

    // Event interface
public:
    virtual void initialize() override final;

    virtual void execute() override final;

    virtual void reset() override final;

};

template<class rT, typename pT>
void KMCSolver<rT, pT>::setSelectedReaction(rT &reaction)
{
    m_selectedReaction = &reaction;
}


template<class rT, typename pT>
KMCSolver<rT, pT>::KMCSolver()
{

}

template<class rT, typename pT>
KMCSolver<rT, pT>::~KMCSolver()
{

}

template<class rT, typename pT>
double KMCSolver<rT, pT>::currentTimeStep() const
{
    return m_currentTimeStep;
}

template<class rT, typename pT>
void KMCSolver<rT, pT>::setCurrentTimeStep(double currentTimeStep)
{
    m_currentTimeStep = currentTimeStep;
}

template<class rT, typename pT>
void KMCSolver<rT, pT>::updateTime()
{
    setCurrentTimeStep(-std::log(rng.uniform())/getTotalRate());

    BADAss(currentTimeStep(), >, 0, "timestep should be posiive.");

    m_currentTime += currentTimeStep();
}

template<class rT, typename pT>
void KMCSolver<rT, pT>::initialize()
{
    initializeReactions();

    m_currentTime = 0;

    updateTime();
}

template<class rT, typename pT>
void KMCSolver<rT, pT>::execute()
{
    selectReaction();
}

template<class rT, typename pT>
void KMCSolver<rT, pT>::reset()
{
    executeReaction(*m_selectedReaction);

    updateReactions();

    updateTime();
}

} //End of namespace kMC
