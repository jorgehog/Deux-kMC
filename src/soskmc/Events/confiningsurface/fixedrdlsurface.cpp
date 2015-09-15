#include "fixedrdlsurface.h"

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
