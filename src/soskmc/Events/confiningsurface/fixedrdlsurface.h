#pragma once

#include "fixedsurface.h"
#include "rdlsurface.h"

class FixedRDLSurface : public FixedSurface, public RDLSurface
{
public:
    FixedRDLSurface(SOSSolver &solver,
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

    bool hasSurface() const
    {
        return true;
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


    // kMC::Observer interface
public:
    void initializeObserver(const Subjects &subject)
    {
        (void) subject;

        RDLSurface::recalculateAllRDLEnergies();
    }

    void notifyObserver(const Subjects &subject);

};
