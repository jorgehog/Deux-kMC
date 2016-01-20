#include "rdlpotential.h"

#include "rdlsurface.h"


RDLPotential::RDLPotential(SOSSolver &solver, RDLSurface &surface) :
    LocalPotential(solver),
    m_surface(surface)
{

}

double RDLPotential::potential(const uint x, const uint y) const
{
    return m_surface.RDLEnergy(x, y);
}
