#pragma once

#include "diffusion.h"

#include <functional>

using std::function;

class ConcentrationProfile : public Diffusion
{
public:

    using funcType = function<double(const uint, const uint)>;

    ConcentrationProfile(SOSSolver &solver, funcType profileFunction);
    virtual ~ConcentrationProfile();

    static double constantDepositionRate()
    {
        return 1.0;
    }

private:

    const funcType m_profileFunction;

    // Event interface
public:
    void execute() {}

    // Diffusion interface
public:
    double depositionRate(const uint x, const uint y) const;
    uint dissolutionPaths(const uint x, const uint y) const;
    void executeDiffusionReaction(SOSDiffusionReaction *reaction, const int x, const int y, const int z);
    void executeConcentrationBoundaryReaction(ConcentrationBoundaryReaction *reaction);
    bool isBlockedPosition(const uint x, const uint y, const int z) const;

    // HeightConnecter interface
public:
    void setupInitialConditions() {}
    void registerHeightChange(const uint x,
                              const uint y,
                              const int value,
                              std::vector<DissolutionDeposition *> &affectedSurfaceReactions,
                              const uint n);
};

