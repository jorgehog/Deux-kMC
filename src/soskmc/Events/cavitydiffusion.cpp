#include "cavitydiffusion.h"

CavityDiffusion::CavityDiffusion(SolidOnSolidSolver &solver,
                                 const double D,
                                 const uint pointsPerLatticeUnit) :
    SolidOnSolidEvent(solver, "Diffusion"),
    m_D(D),
    m_pointsPerLatticeUnit(pointsPerLatticeUnit)
{
    solver.setDiffusionEvent(*this);
}

void CavityDiffusion::execute()
{

}

void CavityDiffusion::registerHeightChange(const uint x, const uint y)
{
    (void) x;
    (void) y;
}

double CavityDiffusion::localSurfaceSupersaturation(const uint x, const uint y)
{
    (void) x;
    (void) y;

    return 1.0;
}

