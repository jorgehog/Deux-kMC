#include "fixedrdlsurface.h"

#include "../../sossolver.h"

FixedRDLSurface::FixedRDLSurface(SOSSolver &solver,
                                 const double s0,
                                 const double ld,
                                 const double height) :
    ConfiningSurface(solver, "fixedRDLSurface"),
    FixedSurface(solver, height),
    RDLSurface(solver, 0, s0, ld)
{

}

FixedRDLSurface::~FixedRDLSurface()
{

}

void FixedRDLSurface::notifyObserver(const Subjects &subject)
{
    const uint &x = solver().currentSurfaceChange().x;
    const uint &y = solver().currentSurfaceChange().y;

    RDLSurface::recalculateRDLEnergy(x, y);

    FixedSurface::notifyObserver(subject);
}
