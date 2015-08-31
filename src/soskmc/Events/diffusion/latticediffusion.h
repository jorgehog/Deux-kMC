#pragma once

#include "diffusion.h"

#include <unordered_map>
using std::unordered_map;

class LatticeDiffusion : public virtual Diffusion
{
public:

    typedef unordered_map<int, SOSDiffusionReaction*> zmap_type;
    typedef unordered_map<uint, zmap_type> zymap_type;
    typedef unordered_map<uint, zymap_type> map_type;

    LatticeDiffusion(SOSSolver &solver);
    ~LatticeDiffusion();

    SOSDiffusionReaction *addDiffusionReactant(const uint x, const uint y, const int z, bool setRate = true);

    void removeDiffusionReactant(SOSDiffusionReaction *reaction, bool _delete = true);

    void removeDiffusionReactant(const uint x, const uint y, const int z, bool _delete = true);

    SOSDiffusionReaction *diffusionReaction(const uint x, const uint y, const int z) const;

    void clearDiffusionReactions();

    uint  numberOfDiffusionReactions() const
    {
        return m_diffusionReactions.size();
    }

    void deleteQueuedReactions();

    void attachToSurface(const uint x, const uint y, const int z, SOSDiffusionReaction *reaction);

    void registerAffectedAround(const uint x, const uint y, const int z);

    void registerAffectedAroundSingle(const int neighbor, const uint xi, const uint dim, const int z);

    void dumpDiffusingParticles(const uint frameNumber) const;

    void moveReaction(SOSDiffusionReaction *reaction, const uint x, const uint y, const int z);

private:

    vector<SOSDiffusionReaction*> m_diffusionReactions;

    map_type m_diffusionReactionsMap;

    vector<SOSDiffusionReaction*> m_deleteQueue;

    SOSSolver &m_mutexSolver;


    // Event interface
public:
    void execute();
    void reset();

    // Diffusion interface
public:
    virtual void dump(const uint frameNumber) const;
    void setupInitialConditions();
    double depositionRate(const uint x, const uint y) const;
    uint dissolutionPaths(const uint x, const uint y) const;
    void registerHeightChange(const uint x, const uint y, const int delta);
    void executeDiffusionReaction(SOSDiffusionReaction *reaction, const int x, const int y, const int z);
    void executeConcentrationBoundaryReaction(ConcentrationBoundaryReaction *reaction);
    bool isBlockedPosition(const uint x, const uint y, const int z) const;

};

