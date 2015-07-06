#include "fixedsurface.h"

FixedSurface::FixedSurface(SolidOnSolidSolver &solver,
                           const double height) :
    ConfiningSurface(solver, "fixedSurface")
{
    setHeight(height);
}

FixedSurface::~FixedSurface()
{

}

