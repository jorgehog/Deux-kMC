#pragma once

#include "offlatticemontecarlo.h"


class FirstPassageContinuum : public OfflatticeMonteCarlo
{
public:
    typedef std::function<double(const FirstPassageContinuum *, uint, uint, uint)> rate_func_type;

    FirstPassageContinuum(SOSSolver &solver,
                          const double maxdt,
                          rate_func_type rateLaw);
    ~FirstPassageContinuum();

private:
    const rate_func_type m_rateLaw;

    // Event interface
public:
    virtual void execute();

    // OfflatticeMonteCarlo interface
public:
    virtual double calculateLocalRate(const uint x, const uint y, const uint n) const
    {
        return m_rateLaw(this, x, y, n);
    }
};

