#pragma once

#include "boundary1d.h"

namespace kMC
{

class FiniteSize : public Boundary1D
{
public:
    FiniteSize(const Boundary::orientations orientation);
    ~FiniteSize();

    // Boundary interface
public:
    double transformCoordinate(const double xi) const;
    void closestImage(const double xi, const double xti, double &dxi) const;
};

}
