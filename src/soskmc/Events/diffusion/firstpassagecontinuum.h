#pragma once

#include "offlatticemontecarlo.h"


class FirstPassageContinuum : public OfflatticeMonteCarlo
{
public:
    FirstPassageContinuum(SOSSolver &solver,
                          const double maxdt,
                          const double n,
                          const double c);

    ~FirstPassageContinuum();

private:
    double m_n;
    double m_c;

    // Event interface
public:
    virtual void execute();

    // OfflatticeMonteCarlo interface
public:
    virtual double calculateLocalRateOverD(const uint x, const uint y, const uint n) const;
};

