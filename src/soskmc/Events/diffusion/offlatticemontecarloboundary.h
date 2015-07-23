#pragma once

#include "offlatticemontecarlo.h"

class SOSDiffusionReaction;

class OfflatticeMonteCarloBoundary : public OfflatticeMonteCarlo
{
public:
    OfflatticeMonteCarloBoundary(SOSSolver &solver,
                                 const double dt,
                                 const uint boundarySpacing = 3);
    ~OfflatticeMonteCarloBoundary();

    bool checkIfEnoughRoom() const;

    void addDiffusionReactant(const uint x, const uint y, const int z);

    void removeDiffusionReactant(SOSDiffusionReaction *reaction);

private:

    const uint m_boundarySpacing;

    vector<SOSDiffusionReaction*> m_diffusionReactions;

    SOSSolver &m_mutexSolver;

    // Event interface
public:
    void execute();
    void initialize();
    void reset();

    // Diffusion interface
public:
    void setupInitialConditions();
    double depositionRate(const uint x, const uint y) const;
    void registerHeightChange(const uint x, const uint y, const int delta);

    // OfflatticeMonteCarlo interface
public:
    void onInsertParticle(const double x, const double y, const double z);
    void onRemoveParticle(const uint n);
};

