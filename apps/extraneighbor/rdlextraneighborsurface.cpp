#include "rdlextraneighborsurface.h"

RDLExtraNeighborSurface::RDLExtraNeighborSurface(SOSSolver &solver,
                                                 RDLPotential &rdlPotential,
                                                 ExtraNeighbor &extraNeighborPotential,
                                                 const double Pl) :
    ConfiningSurface(solver, "RDLExtraNeighborSurface"),
    m_rdlPotential(rdlPotential),
    m_extraNeighborPotential(extraNeighborPotential),
    m_Pl(Pl)
{
    registerObserver(&rdlPotential);
    registerObserver(&extraNeighborPotential);
}

void RDLExtraNeighborSurface::findNewHeight()
{
    const double eps = 1E-10;

    double h = height();
    double hPrev;

    cout << h << endl;

    do
    {
        hPrev = h;

        h = iteratingExpression(hPrev);

    } while (abs(h - hPrev) > eps);

    cout << h << endl;

    setHeight(h);
}

double RDLExtraNeighborSurface::iteratingExpression(const double hl) const
{
    double nFactor = 6*ExtraNeighbor::m_scaling*m_rdlPotential.lD()*(1./solver().area() + m_rdlPotential.s0()/(1 - exp(-1/m_rdlPotential.lD())));

    double nSum = 0;
    double theta = 0;

    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            const int &hi = solver().height(x, y);

            theta += exp(hi/m_rdlPotential.lD());

            if (hl - hi < 2 && hl > hi)
            {
                nSum += hl - hi;
            }

        }
    }

    theta /= solver().area();

    return 1 - m_rdlPotential.lD()*log((m_Pl+nSum*nFactor)/(m_rdlPotential.s0()*theta));
}

void RDLExtraNeighborSurface::initializeObserver(const Subjects &subject)
{
    (void) subject;

    setHeight(solver().heights().max() + 1);
    findNewHeight();
}

void RDLExtraNeighborSurface::notifyObserver(const Subjects &subject)
{
    (void) subject;

    findNewHeight();
}

void RDLExtraNeighborSurface::execute()
{

}

bool RDLExtraNeighborSurface::hasSurface() const
{
    return true;
}

bool RDLExtraNeighborSurface::acceptDiffusionMove(const double x0, const double y0, const double z0, const double x1, const double y1, const double z1) const
{
    (void) x0;
    (void) y0;
    (void) z0;
    (void) x1;
    (void) y1;
    (void) z1;

    return true;
}

double RDLExtraNeighborSurface::diffusionDrift(const double x, const double y, const double z) const
{
    (void) x;
    (void) y;
    (void) z;

    return 0.0;
}
