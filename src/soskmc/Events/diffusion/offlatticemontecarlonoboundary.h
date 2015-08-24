#pragma once

#include "offlatticemontecarlo.h"

class OfflatticeMonteCarloNoBoundary : public OfflatticeMonteCarlo
{
public:
    OfflatticeMonteCarloNoBoundary(SOSSolver &solver,
                                   const double dt);

    ~OfflatticeMonteCarloNoBoundary();

    double calculateTimeStep(const double initialCondition, bool calculateDissolutionRate = false);

    double depositionRate(const uint x, const uint y, double timeStep) const;

    double calculateLocalProbability(const uint x, const uint y, const uint n, const double timeStep) const;

    double calculateLocalProbability(const uint x, const uint y, const uint n) const
    {
        return calculateLocalProbability(x, y, n, m_currentTimeStep);
    }


private:

    double m_currentTimeStep;

    static constexpr double m_eps = 0.00001;

    cube m_localProbabilities;



    // Event interface
public:
    void execute();

    void reset();

    // Diffusion interface
public:
    void setupInitialConditions();
    double depositionRate(const uint x, const uint y) const
    {
        return depositionRate(x, y, m_currentTimeStep);
    }
    void registerHeightChange(const uint x, const uint y, const int delta);
    void executeDiffusionReaction(SOSDiffusionReaction *reaction, const int x, const int y, const int z);
    void executeConcentrationBoundaryReaction(ConcentrationBoundaryReaction *reaction);
    bool isBlockedPosition(const uint x, const uint y, const int z) const;
    void dump(const uint frameNumber) const;

    // OfflatticeMonteCarlo interface
public:
    void onInsertParticle(const double x, const double y, const double z);
    void onRemoveParticle(const uint n);

};
