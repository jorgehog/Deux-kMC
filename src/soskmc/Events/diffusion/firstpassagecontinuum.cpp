#include "firstpassagecontinuum.h"

#include "../../sossolver.h"

FirstPassageContinuum::FirstPassageContinuum(SOSSolver &solver,
                                             const double maxdt, const double c) :
    Diffusion(solver, "FirstPassage", "", true),
    OfflatticeMonteCarlo(solver, maxdt),
    m_c(c)
{

}

FirstPassageContinuum::~FirstPassageContinuum()
{

}

void FirstPassageContinuum::setc(const double c)
{
    m_c = c;

    calculateLocalRates();
}

void FirstPassageContinuum::execute()
{
    setValue(acceptanceRatio());
}

double FirstPassageContinuum::calculateLocalRateOverD(const uint x, const uint y, const uint n) const
{
    const int z = solver().height(x, y) + 1;

    const double dr2 = solver().closestSquareDistance(x, y, z,
                                                      particlePositions(0, n),
                                                      particlePositions(1, n),
                                                      particlePositions(2, n));

    return calculateLocalRateOverD(dr2);
}

double FirstPassageContinuum::calculateLocalRateOverD(const double rSquared) const
{
    return m_c/(rSquared*rSquared);
}
