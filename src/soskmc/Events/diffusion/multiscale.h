#pragma once

#include "offlatticemontecarlo.h"
#include "latticediffusion.h"

class SOSDiffusionReaction;

class OfflatticeMonteCarloBoundary : public OfflatticeMonteCarlo, public LatticeDiffusion
{
public:
    OfflatticeMonteCarloBoundary(SOSSolver &solver,
                                 const double maxdt,
                                 const uint boundarySpacing = 3);
    ~OfflatticeMonteCarloBoundary();

    bool checkIfEnoughRoom() const;


private:

    const uint m_boundarySpacing;

    // Event interface
public:
    void execute();
    void reset()
    {
        OfflatticeMonteCarlo::reset();
    }

    // Diffusion interface
public:

    double depositionRate(const uint x, const uint y) const
    {
        return LatticeDiffusion::depositionRate(x, y);
    }

    uint dissolutionPaths(const uint x, const uint y) const
    {
        return OfflatticeMonteCarlo::dissolutionPaths(x, y);
    }

    void executeDiffusionReaction(SOSDiffusionReaction *reaction, const int x, const int y, const int z);

    void executeConcentrationBoundaryReaction(ConcentrationBoundaryReaction *reaction);

    bool isBlockedPosition(const uint x, const uint y, const int z) const
    {
        return LatticeDiffusion::isBlockedPosition(x, y, z);
    }

    void dump(const uint frameNumber, const string path = "/tmp") const;


    // HeightConnecter interface
public:
    void setupInitialConditions();
    void registerHeightChange(const uint x,
                              const uint y,
                              const int value,
                              std::vector<DissolutionDeposition *> &affectedSurfaceReactions,
                              const uint n)
    {
        LatticeDiffusion::registerHeightChange(x, y, value, affectedSurfaceReactions, n);
    }


    // OfflatticeMonteCarlo interface
public:
    double calculateLocalRate(const uint x, const uint y, const uint n) const;
};

