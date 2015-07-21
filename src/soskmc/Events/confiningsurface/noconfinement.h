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
    void setupInitialConditions();
    void registerHeightChange(const uint x, const uint y, std::vector<DissolutionDeposition *> affectedReactions, const uint n);
    double confinementEnergy(const uint x, const uint y);
    bool acceptDiffusionMove(const double x0, const double y0, const double z0, const double x1, const double y1, const double z1) const;
    double diffusionDrift(const double x, const double y, const double z) const;
};
