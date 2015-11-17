#pragma once

#include "confiningsurface.h"

class FixedSurface : public virtual ConfiningSurface
{
public:
    FixedSurface(SOSSolver &solver,
                 const double height);

    virtual ~FixedSurface();

    // Event interface
public:
    void execute()
    {
        setValue(height());
    }

    // ConfiningSurface interface
public:

    bool hasSurface() const
    {
        return true;
    }

    double confinementEnergy(const uint x, const uint y)
    {
        (void) x;
        (void) y;

        return 0;
    }

    bool acceptDiffusionMove(const double x0, const double y0, const double z0,
                             const double x1, const double y1, const double z1) const;

    double diffusionDrift(const double x, const double y, const double z) const;

    // HeightObserver interface
public:
    void initializeObserver()
    {
        //pass
    }

    virtual void notifyObserver()
    {
        //pass
    }
};
