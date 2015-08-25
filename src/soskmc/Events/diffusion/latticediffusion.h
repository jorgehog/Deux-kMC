#pragma once

#include "diffusion.h"

class LatticeDiffusion : public virtual Diffusion
{
public:
    LatticeDiffusion(SOSSolver &solver);
    ~LatticeDiffusion();

    SOSDiffusionReaction *addDiffusionReactant(const uint x, const uint y, const int z, bool setRate = true);

    void removeDiffusionReactant(SOSDiffusionReaction *reaction, bool _delete = true);

    void removeDiffusionReactant(const uint x, const uint y, const int z, bool _delete = true);

    SOSDiffusionReaction *diffusionReaction(const uint x, const uint y, const int z) const;

    SOSDiffusionReaction *diffusionReaction(const uint n) const;

    void clearDiffusionReactions();

    uint  numberOfDiffusionReactions() const
    {
        return m_diffusionReactions.size();
    }

    void deleteQueuedReactions();

    void attachToSurface(const uint x, const uint y, const int z, SOSDiffusionReaction *reaction);

private:

    vector<SOSDiffusionReaction*> m_diffusionReactions;

    vector<SOSDiffusionReaction*> m_deleteQueue;

    SOSSolver &m_mutexSolver;


    // Event interface
public:
    void execute();


    // Diffusion interface
public:
    virtual void dump(const uint frameNumber) const;
    void setupInitialConditions();
    double depositionRate(const uint x, const uint y) const;
    void registerHeightChange(const uint x, const uint y, const int delta);
    void executeDiffusionReaction(SOSDiffusionReaction *reaction, const int x, const int y, const int z);
    void executeConcentrationBoundaryReaction(ConcentrationBoundaryReaction *reaction);
    bool isBlockedPosition(const uint x, const uint y, const int z) const;

};

