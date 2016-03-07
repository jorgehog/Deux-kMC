#include "firstpassagecontinuum.h"

#include "../../sossolver.h"

FirstPassageContinuum::FirstPassageContinuum(SOSSolver &solver,
                                             const double maxdt,
                                             const int depositionBoxHalfSize,
                                             const double c) :
    Diffusion(solver, "FirstPassage"),
    OfflatticeMonteCarlo(solver, maxdt),
    m_c(c),
    m_depositionBoxHalfSize(depositionBoxHalfSize),
    m_allSimBoxes(solver.length(), solver.width())
{
    int xTrans;
    int yTrans;

    const int &l = depositionBoxHalfSize;

    for (uint x = 0; x < solver.length(); ++x)
    {
        for (uint y = 0; y < solver.width(); ++y)
        {
            m_allSimBoxes(x, y).set_size(2*l + 1, 2*l + 1);

            for (int xscan = -l; xscan <= l; ++xscan)
            {
                for (int yscan = -l; yscan <= l; ++yscan)
                {
                    solver.boundaryLatticeTransform(xTrans, yTrans, (int)x + xscan, (int)y + yscan, 0);

                    m_allSimBoxes(x, y)(xscan + l, yscan + l) = make_pair(xTrans, yTrans);
                }
            }
        }
    }
}

FirstPassageContinuum::~FirstPassageContinuum()
{

}

void FirstPassageContinuum::setc(const double c)
{
    m_c = c;

    calculateLocalRatesAndUpdateDepositionRates();
}


double FirstPassageContinuum::localRateOverD(const double rSquared) const
{
    return m_c/(rSquared*rSquared);
}

void FirstPassageContinuum::checkTrans(int xTrans,
                                       int yTrans,
                                       const int x,
                                       const int y,
                                       const int dx,
                                       const int dy)
{
    int xTrans1;
    int yTrans1;

    solver().boundaryLatticeTransform(xTrans1, yTrans1, x + dx, y + dy, 0);

    if (xTrans != xTrans1 || yTrans != yTrans1)
    {
        cout << "error " << xTrans << " " << xTrans1 << " " << yTrans << " " << yTrans1 << endl;
        exit(1);
    }
}
