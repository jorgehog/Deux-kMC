#include "noconfinement.h"

NoConfinement::NoConfinement(SOSSolver &solver) :
    ConfiningSurface(solver, "NoConfinement")
{

}

NoConfinement::~NoConfinement()
{

}


void NoConfinement::execute()
{

}

void NoConfinement::setupInitialConditions()
{

}

void NoConfinement::registerHeightChange(const uint x, const uint y, std::vector<DiffusionDeposition *> affectedReactions, const uint n)
{
    (void) x;
    (void) y;
    (void) affectedReactions;
    (void) n;
}

double NoConfinement::confinementEnergy(const uint x, const uint y)
{
    (void) x;
    (void) y;

    return 0;
}

bool NoConfinement::acceptDiffusionMove(const double x0, const double y0, const double z0,
                                        const double x1, const double y1, const double z1) const
{
    (void) x0;
    (void) y0;
    (void) z0;
    (void) x1;
    (void) y1;
    (void) z1;

    BADAssBreak("Diffusion is not supported without confinement.");
    return 0;
}

double NoConfinement::diffusionDrift(const double x, const double y, const double z) const
{
    (void) x;
    (void) y;
    (void) z;

    BADAssBreak("Diffusion is not supported without confinement.");
    return 0;
}
