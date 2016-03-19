#include "radialfirstpassage.h"

#include "../../sossolver.h"
#include "dissolutiondeposition.h"
#include "../confiningsurface/confiningsurface.h"

RadialFirstPassage::RadialFirstPassage(SOSSolver &solver,
                                       const double maxdt,
                                       const int depositionBoxHalfSize,
                                       const double c) :
    Diffusion(solver, "RadialFirstPassage"),
    FirstPassageContinuum(solver, maxdt, depositionBoxHalfSize, c)
{

}

double RadialFirstPassage::localRateOverD(const uint x, const uint y, const uint n) const
{
    const int z = solver().height(x, y) + 1;

    const double dr2 = solver().closestSquareDistance(x, y, z,
                                                      particlePositions(0, n),
                                                      particlePositions(1, n),
                                                      particlePositions(2, n));

    return FirstPassageContinuum::localRateOverD(dr2);
}

void RadialFirstPassage::execute()
{
}

void RadialFirstPassage::calculateLocalRatesAndUpdateDepositionRates()
{
    resetDepositionRates();

    const int &l = depositionBoxHalfSize();

    int xTrans;
    int yTrans;

    SurfaceReaction *r;
    double r2;
    double localRate;

    const int hmax = solver().heights().max();
    const double hUpper = solver().confiningSurfaceEvent().height() - 2;

    const double D = DScaled();

    m_localRates.zeros();
    for (uint n = 0; n < nOfflatticeParticles(); ++n)
    {
        const double &zp = particlePositions(2, n);

        if (zp > hmax + l + 1)
        {
            continue;
        }

        const double &xp = particlePositions(0, n);
        const double &yp = particlePositions(1, n);

        const int ix = roundfast(xp);
        const int iy = roundfast(yp);

        const double dx = ix - xp;
        const double dy = iy - yp;

        for (int xscan = -l; xscan <= l; ++xscan)
        {
            for (int yscan = -l; yscan <= l; ++yscan)
            {
                getTrans(xTrans, yTrans, ix, iy, xscan, yscan);

                const int &h = solver().height(xTrans, yTrans);

                //if there is no room to deposit.
                if (h > hUpper)
                {
                    continue;
                }

                const double dz2 = (h+1-zp)*(h+1-zp);

                if (dz2 < l*l)
                {
                    r2 = dz2 + (xscan + dx)*(xscan + dx) + (yscan + dy)*(yscan+dy);

                    BADAssClose(r2, solver().closestSquareDistance(xTrans, yTrans, h+1, xp, yp, zp), 1E-3);

                    r = &solver().surfaceReaction(xTrans, yTrans);
                    localRate = c()/(r2*r2);
                    m_localRates(xTrans, yTrans, n) = localRate;
                    r->setDepositionRate(r->depositionRate() + localRate*D);
                }
            }
        }
    }
}
