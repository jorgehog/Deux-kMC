#pragma once

#include "boundary.h"

namespace kMC
{

class FiniteSize : public Boundary
{
public:
    FiniteSize(const Boundary::orientations orientation);
    ~FiniteSize();

    // Boundary interface
public:
    double transformCoordinate(const double xi) const;
    std::vector<double> imagesOf(const double xi) const;
};

}
