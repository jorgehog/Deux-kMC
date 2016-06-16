#pragma once

#include "offlatticemontecarlo.h"


inline int roundfast(const double x)
{
    const double xshift = x + 0.5;
    int rval;

    if (xshift < 0)
    {
        rval = int(xshift) - 1;
    }

    else
    {
        rval = int(xshift);
    }

    BADAss(rval, ==, int(round(x)));

    return rval;
}

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

    void getTrans(int &xTrans, int &yTrans,
                  const uint x, const uint y,
                  const int dx, const int dy);

    void checkTrans(int xTrans, int yTrans,
                    const int x, const int y,
                    const int dx, const int dy);

private:
    double m_c;

    const int m_depositionBoxHalfSize;

    field<field<pair<int, int>>> m_allSimBoxes;


};

