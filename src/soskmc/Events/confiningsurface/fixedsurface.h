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
        //pass
    }

    // ConfiningSurface interface
public:
    void setupInitialConditions()
    {
        //pass
    }

    void registerHeightChange(const uint x, const uint y, std::vector<DiffusionDeposition *> affectedReactions, const uint n)
    {
        (void) affectedReactions;
        (void) n;

        if (solver().height(x, y) > height())
        {
            terminateLoop("CRASH CRASH...");
        }
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
};
