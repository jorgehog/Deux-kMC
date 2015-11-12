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

private:

    const funcType m_profileFunction;

    // Event interface
public:
    virtual void execute() {}

    // Diffusion interface
public:
    virtual double depositionRate(const uint x, const uint y) const;
    virtual uint dissolutionPaths(const uint x, const uint y) const;
    virtual void executeDiffusionReaction(SOSDiffusionReaction *reaction, const int x, const int y, const int z);
    virtual void executeConcentrationBoundaryReaction(ConcentrationBoundaryReaction *reaction);
    virtual bool isBlockedPosition(const uint x, const uint y, const int z) const;
    double concentration() const;

    // HeightConnecter interface
public:
    virtual void setupInitialConditions() {}
    virtual void registerHeightChange(const uint x,
                              const uint y,
                              const int value,
                              std::vector<DissolutionDeposition *> &affectedSurfaceReactions,
                              const uint n);
};

