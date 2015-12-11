#include "radialfirstpassage.h"

#include "../../sossolver.h"
#include "dissolutiondeposition.h"


RadialFirstPassage::RadialFirstPassage(SOSSolver &solver,
                                       const double maxdt,
                                       const int depositionBoxHalfSize,
                                       const double c) :
    Diffusion(solver, "AStarFirstPassage"),
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

void RadialFirstPassage::calculateLocalRates()
{
    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            solver().surfaceReaction(x, y).setDepositionRate(0);
        }
    }

    const int &l = depositionBoxHalfSize();

    int xTrans;
    int yTrans;

    DissolutionDeposition *r;
    double r2;
    double localRate;

    const int hmax = solver().heights().max();

    const double D = DScaled();

    for (uint n = 0; n < nOfflatticeParticles(); ++n)
    {

        for (uint x = 0; x < solver().length(); ++x)
        {
            for (uint y = 0; y < solver().width(); ++y)
            {
                m_localRates(x, y, n) = 0;
            }
        }

        const double &zp = particlePositions(2, n);

        if (zp > hmax + l + 1)
        {
            continue;
        }

        const double &xp = particlePositions(0, n);
        const double &yp = particlePositions(1, n);

        const int ix = int(round(xp));
        const int iy = int(round(yp));
        const int iz = int(round(zp));

        const double dx = ix - xp;
        const double dy = iy - yp;

        for (int xscan = -l; xscan <= l; ++xscan)
        {
            for (int yscan = -l; yscan <= l; ++yscan)
            {
                solver().boundaryLatticeTransform(xTrans, yTrans, ix + xscan, iy + yscan, iz);

                const int &h = solver().height(xTrans, yTrans);

                const double dz2 = (h+1-zp)*(h+1-zp);

                if (dz2 < l*l)
                {
                    r2 = dz2 + (xscan + dx)*(xscan + dx) + (yscan + dy)*(yscan+dy);

                    BADAssClose(r2, solver().closestSquareDistance(xTrans, yTrans, h+1, xp, yp, zp), 1E-3);

                    r = &solver().surfaceReaction(xTrans, yTrans);
                    localRate = FirstPassageContinuum::localRateOverD(r2);
                    m_localRates(xTrans, yTrans, n) = localRate;
                    r->setDepositionRate(r->depositionRate() + localRate*D);
                }
            }
        }
    }
}
