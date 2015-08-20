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
    int transformCoordinate(const int xi) const;
};

}
