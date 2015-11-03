#pragma once

#include "kmcassets.h"

#include "RNG/rng.h"

#include <ignis/include/ignis.h>

#include <unordered_set>

using std::unordered_set;

namespace kMC
{

class Reaction;
class Boundary;

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

    void setBoundaries(vector<vector<const Boundary *> > boundaries)
    {
        m_boundaries = boundaries;
    }

    const Boundary *boundary(const uint i, const uint j) const
    {
        return m_boundaries.at(i).at(j);
    }

    void initializeReactions();

    void updateAffectedReactions();

    Reaction *getReaction(const uint n) const
    {
        BADAss(n, <, m_reactions.size(), "index out of bounds.");

        return m_reactions[n];
    }

    const double &nextRandomLogNumber() const
    {
        return m_nextRandomLogNumber;
    }

    void addReaction(Reaction *reaction);

    void removeReaction(const uint n);

    void removeReaction(Reaction *reaction);

    void removeIfAffected(Reaction *reaction);

    uint numberOfReactions() const
    {
        return m_reactions.size();
    }

    void registerAffectedReaction(Reaction *reaction)
    {
        BADAss(std::find(m_reactions.begin(), m_reactions.end(), reaction), !=, m_reactions.end(), "affecting reaction which is not in reaction list.");
        m_affectedReactions.insert(reaction);
    }

    const unordered_set<Reaction*>  &affectedReactions() const
    {
        return m_affectedReactions;
    }

private:
    vector<Reaction*> m_reactions;

    unordered_set<Reaction*> m_affectedReactions;

    vector<vector<const Boundary*>> m_boundaries;

    Reaction *m_selectedReaction;

    double m_currentTimeStep;
    double m_currentTime;

    double m_totalRate;

    double m_nextRandomLogNumber;

    vector<double> m_cumsumRates;

    double getRandomLogNumber() const;

    void updateTime();

    void setCurrentTimeStep(double currentTimeStep);

    void resizeCumsumRates();

    void getCumsumAndTotalRate();

    // Event interface
public:
    virtual void initialize();

    virtual void execute();

    virtual void reset();

};

} //End of namespace kMC
