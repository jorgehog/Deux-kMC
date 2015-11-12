#pragma once

#include "fixedsurface.h"

class ConstantConfinement : public FixedSurface
{
public:
    ConstantConfinement(SOSSolver &solver, const double height);
    ~ConstantConfinement();

private:

    double m_h;

    // HeightConnecter interface
public:
    void registerHeightChange(const uint x, const uint y, const int value, std::vector<DissolutionDeposition *> &affectedSurfaceReactions, const uint nAffectedSurfaceReactions);

    // Event interface
public:
    void initialize();
};

