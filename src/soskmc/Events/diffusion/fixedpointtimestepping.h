#pragma once

#include "offlatticemontecarlo.h"

class FixedPointTimeStepping : public OfflatticeMonteCarlo
{
public:
    FixedPointTimeStepping(SOSSolver &solver, const double maxdt);

    double calculateTimeStep(const double initialCondition, bool calculateDissolutionRate = false);

    double calculateLocalRate(const uint x, const uint y, const uint n, const double timeStep) const;


private:
    SOSSolver &m_mutexSolver;

    double m_currentTimeStep;

    static constexpr double m_eps = 0.00001;


    // Event interface
public:
    void execute();


    // Diffusion interface
public:
    void setupInitialConditions();
    void registerHeightChange(const uint x, const uint y, const int delta);
    double calculateLocalRate(const uint x, const uint y, const uint n) const;

};

