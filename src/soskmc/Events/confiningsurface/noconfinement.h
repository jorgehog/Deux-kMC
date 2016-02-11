#pragma once

#include "confiningsurface.h"

class NoConfinement : public ConfiningSurface
{
public:
    NoConfinement(SOSSolver &solver);
    ~NoConfinement();

    // Event interface
public:
    void execute();

    // ConfiningSurface interface
public:
    bool hasSurface() const
    {
        return false;
    }
    bool acceptDiffusionMove(const double x0,
                             const double y0,
                             const double z0,
                             const double x1,
                             const double y1,
                             const double z1) const;
    double diffusionDrift(const double x,
                          const double y,
                          const double z) const;

    // kMC::Observer interface
public:
    void initializeObserver(const Subjects &subject);
    void notifyObserver(const Subjects &subject);
};
