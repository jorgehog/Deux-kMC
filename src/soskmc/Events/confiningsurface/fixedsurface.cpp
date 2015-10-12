#include "fixedsurface.h"

#include "sossolver.h"

FixedSurface::FixedSurface(SOSSolver &solver,
                           const double height) :
    ConfiningSurface(solver, "fixedSurface")
{
    setHeight(height);
}

FixedSurface::~FixedSurface()
{

}

bool FixedSurface::acceptDiffusionMove(const double x0, const double y0, const double z0, const double x1, const double y1, const double z1) const
{
    (void) x0;
    (void) y0;
    (void) z0;
    (void) x1;
    (void) y1;
    (void) z1;

    return true;
}

double FixedSurface::diffusionDrift(const double x, const double y, const double z) const
{
    (void) x;
    (void) y;
    (void) z;

    return 0;
}

void FixedSurface::registerHeightChange(const uint x, const uint y, const int value, std::vector<DissolutionDeposition *> &affectedSurfaceReactions, const uint n)
{
    (void) value;
    (void) affectedSurfaceReactions;
    (void) n;

    if (solver().height(x, y) > height())
    {
        terminateLoop("CRASH CRASH...");
    }
}
