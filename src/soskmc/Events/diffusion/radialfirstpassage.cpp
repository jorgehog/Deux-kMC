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

double RadialFirstPassage::_localRateOverD(const uint x, const uint y, const uint n) const
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
    const int &l = depositionBoxHalfSize();

    int xTrans;
    int yTrans;

    double r2;

    const int hmax = solver().heights().max();
    const double hUpper = solver().confiningSurfaceEvent().height() - 2;

    const double D = DScaled();

    for (uint n = 0; n < nOfflatticeParticles(); ++n)
    {
        resetLocalRates(n);

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

                if (solver().isOutsideBox(xTrans, yTrans))
                {
                    continue;
                }

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

                    //for boundaries such as reflective we cannot use shortest image.
                    //BADAssClose(r2, solver().closestSquareDistance(xTrans, yTrans, h+1, xp, yp, zp), 1E-3);
                    m_localRates(xTrans, yTrans, n) += localRateOverD(r2);
                }
            }
        }
    }

    for (SurfaceReaction *r : solver().surfaceReactions())
    {
        double localrate = 0;
        for (uint n = 0; n < nOfflatticeParticles(); ++n)
        {
            localrate += m_localRates(r->x(), r->y(), n);
        }

        r->setDepositionRate(localrate*D);
    }
}
