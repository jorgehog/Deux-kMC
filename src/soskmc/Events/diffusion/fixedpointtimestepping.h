#pragma once

#include "offlatticemontecarlo.h"

class FixedPointTimeStepping : public OfflatticeMonteCarlo
{
public:
    FixedPointTimeStepping(SOSSolver &solver, const double maxdt);

    void calculateTimeStep(bool calculateDissolutionRate = false);

    double calculateLocalRateOverD(const uint x, const uint y, const uint n, const double timeStep) const;

    double calculateLocalRateOverD(const double rSquared) const {(void)rSquared;return 0;}

private:
    SOSSolver &m_mutexSolver;

    double m_currentTimeStep;

    static constexpr double m_eps = 0.00001;


    // Event interface
public:
    void execute();


    // Diffusion interface
public:
    double calculateLocalRateOverD(const uint x, const uint y, const uint n) const;


    // kMC::Observer interface
public:
    void initializeObserver(const Subjects &subject);
    void notifyObserver(const Subjects &subject);

    // OfflatticeMonteCarlo interface
public:
    void calculateLocalRates();
};

