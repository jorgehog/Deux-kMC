#pragma once

#include "offlatticemontecarlo.h"


class FirstPassageContinuum : public OfflatticeMonteCarlo
{
public:
    FirstPassageContinuum(SOSSolver &solver,
                          const double maxdt,
                          const double c);

    ~FirstPassageContinuum();

    void setc(const double c);

    const double &c() const
    {
        return m_c;
    }

private:
    double m_c;

    // Event interface
public:
    virtual void execute();

    // OfflatticeMonteCarlo interface
public:
    virtual double calculateLocalRateOverD(const uint x, const uint y, const uint n) const;
    virtual double calculateLocalRateOverD(const double rSquared) const;
};

