#pragma once

#include "firstpassagecontinuum.h"

class RadialFirstPassage : public FirstPassageContinuum
{
public:
    RadialFirstPassage(SOSSolver &solver,
                       const double maxdt,
                       const int depositionBoxHalfSize,
                       const double c);

    //not used
    double localRateOverD(const uint x, const uint y, const uint n) const;

    // Event interface
public:
    void execute();

    // OfflatticeMonteCarlo interface
public:
    void calculateLocalRatesAndUpdateDepositionRates();
};

