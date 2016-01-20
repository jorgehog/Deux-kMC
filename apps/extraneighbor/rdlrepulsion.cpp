#include "rdlrepulsion.h"

RDLRepulsion::RDLRepulsion(SOSSolver &solver, const double s0, const double ld, const double Pl) :
    LocalCachedPotential(solver),
    m_s0(s0),
    m_ld(ld),
    m_Pl(Pl)
{

}

void RDLRepulsion::notifyObserver(const Subjects &subject)
{

}

double RDLRepulsion::potentialFunction(const uint x, const uint y) const
{
    const double &hl = solver().confiningSurfaceEvent().height();
    const int &hi = solver().height(x, y);

    const double dh = hl - hi;

    return m_s0*exp(-dh/m_ld);
}


