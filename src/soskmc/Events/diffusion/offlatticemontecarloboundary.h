#pragma once

#include "offlatticemontecarlo.h"

class OfflatticeMonteCarloBoundary : public OfflatticeMonteCarlo
{
public:
    OfflatticeMonteCarloBoundary(SolidOnSolidSolver &solver,
                                 const double dt,
                                 const uint boundarySpacing = 3);
    ~OfflatticeMonteCarloBoundary();

    bool checkIfEnoughRoom() const;

private:

    const uint m_boundarySpacing;

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

