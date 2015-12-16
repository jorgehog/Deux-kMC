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

    if (solver().currentSurfaceChange().type == ChangeTypes::Double)
    {
        const uint &x1 = solver().currentSurfaceChange().x1;
        const uint &y1 = solver().currentSurfaceChange().y1;

        RDLSurface::recalculateRDLEnergy(x1, y1);
    }

    FixedSurface::notifyObserver(subject);
}
