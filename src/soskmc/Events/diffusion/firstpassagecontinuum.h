#pragma once

#include "offlatticemontecarlo.h"


class FirstPassageContinuum : public OfflatticeMonteCarlo
{
public:
    FirstPassageContinuum(SOSSolver &solver,
                          const double maxdt,
                          const int depositionBoxHalfSize,
                          const double c);

    ~FirstPassageContinuum();

    void setc(const double c);

    const double &c() const
    {
        return m_c;
    }

    const int &depositionBoxHalfSize() const
    {
        return m_depositionBoxHalfSize;
    }

    double localRateOverD(const double rSquared) const;

private:
    double m_c;

    const int m_depositionBoxHalfSize;

};

