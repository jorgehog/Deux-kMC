#pragma once

#include "fixedsurface.h"
#include "rdlsurface.h"

class FixedRDLSurface : public FixedSurface, public RDLSurface
{
public:
    FixedRDLSurface(SOSSolver &solver,
                    const double E0,
                    const double s0,
                    const double ld,
                    const double height);

    ~FixedRDLSurface();

    // Event interface
public:
    void execute()
    {
        return FixedSurface::execute();
    }
    void reset()
    {
        return RDLSurface::_validateStoredEnergies();
    }

    // ConfiningSurface interface
public:
    void setupInitialConditions()
    {
        RDLSurface::recalculateAllRDLEnergies();
    }

    void registerHeightChange(const uint x, const uint y, std::vector<DiffusionDeposition *> affectedReactions, const uint n)
    {
        (void) affectedReactions;
        (void) n;

        RDLSurface::recalculateRDLEnergy(x, y);

        FixedSurface::registerHeightChange(x, y, affectedReactions, n);
    }

    double confinementEnergy(const uint x, const uint y)
    {
        return RDLSurface::confinementEnergy(x, y);
    }

    double diffusionDrift(const double x, const double y, const double z) const
    {
        return RDLSurface::diffusionDrift(x, y, z);
    }

    bool acceptDiffusionMove(const double x0, const double y0, const double z0,
                             const double x1, const double y1, const double z1) const
    {
        return RDLSurface::acceptDiffusionMove(x0, y0, z0, x1, y1, z1);
    }
};
