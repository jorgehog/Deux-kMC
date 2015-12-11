#include "firstpassagecontinuum.h"

#include "../../sossolver.h"

FirstPassageContinuum::FirstPassageContinuum(SOSSolver &solver,
                                             const double maxdt,
                                             const int depositionBoxHalfSize,
                                             const double c) :
    Diffusion(solver, "FirstPassage"),
    OfflatticeMonteCarlo(solver, maxdt),
    m_c(c),
    m_depositionBoxHalfSize(depositionBoxHalfSize)
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


double FirstPassageContinuum::localRateOverD(const double rSquared) const
{
    return m_c/(rSquared*rSquared);
}
