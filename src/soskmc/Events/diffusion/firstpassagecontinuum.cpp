#include "firstpassagecontinuum.h"

#include "../../sossolver.h"

FirstPassageContinuum::FirstPassageContinuum(SOSSolver &solver,
                                             const double maxdt,
                                             const double n, const double c) :
    Diffusion(solver, "FirstPassage", "", true),
    OfflatticeMonteCarlo(solver, maxdt),
    m_n(n),
    m_c(c)
{

}

FirstPassageContinuum::~FirstPassageContinuum()
{

}

void FirstPassageContinuum::execute()
{
    setValue(acceptanceRatio());
}

double FirstPassageContinuum::calculateLocalRate(const uint x, const uint y, const uint n) const
{
    const int z = solver().height(x, y) + 1;

    const double dr2 = solver().closestSquareDistance(x, y, z,
                                                      particlePositions(0, n),
                                                      particlePositions(1, n),
                                                      particlePositions(2, n));

    return m_c/std::pow(dr2, m_n/2);
}
