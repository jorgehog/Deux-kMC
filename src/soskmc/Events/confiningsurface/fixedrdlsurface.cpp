#include "fixedrdlsurface.h"

FixedRDLSurface::FixedRDLSurface(SolidOnSolidSolver &solver,
                                 const double E0,
                                 const double s0,
                                 const double ld,
                                 const double height) :
    ConfiningSurface(solver, "fixedRDLSurface"),
    FixedSurface(solver, height),
    RDLSurface(solver, E0, s0, ld)
{

}

FixedRDLSurface::~FixedRDLSurface()
{

}

