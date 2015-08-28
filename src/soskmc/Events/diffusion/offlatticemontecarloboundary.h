#pragma once

#include "offlatticemontecarlo.h"
#include "latticediffusion.h"

class SOSDiffusionReaction;

class OfflatticeMonteCarloBoundary : public OfflatticeMonteCarlo, public LatticeDiffusion
{
public:
    OfflatticeMonteCarloBoundary(SOSSolver &solver,
                                 const double dt,
                                 const uint boundarySpacing = 3);
    ~OfflatticeMonteCarloBoundary();

    bool checkIfEnoughRoom() const;


private:

    const uint m_boundarySpacing;

    // Event interface
public:
    void execute();

    // Diffusion interface
public:
    void setupInitialConditions();

    double depositionRate(const uint x, const uint y) const
    {
        return LatticeDiffusion::depositionRate(x, y);
    }

    uint dissolutionPaths(const uint x, const uint y) const
    {
        return OfflatticeMonteCarlo::dissolutionPaths(x, y);
    }

    void registerHeightChange(const uint x, const uint y, const int delta)
    {
        LatticeDiffusion::registerHeightChange(x, y, delta);
    }

    void executeDiffusionReaction(SOSDiffusionReaction *reaction, const int x, const int y, const int z)
    {
        LatticeDiffusion::executeDiffusionReaction(reaction, x, y, z);
    }

    void executeConcentrationBoundaryReaction(ConcentrationBoundaryReaction *reaction);

    bool isBlockedPosition(const uint x, const uint y, const int z) const
    {
        return LatticeDiffusion::isBlockedPosition(x, y, z);
    }

    void dump(const uint frameNumber) const;

    // OfflatticeMonteCarlo interface
public:
    void onInsertParticle(const double x, const double y, const double z);
    void onRemoveParticle(const uint n);

};

