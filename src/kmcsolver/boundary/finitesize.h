#pragma once

#include "boundary.h"

namespace kMC
{

class FiniteSize : public Boundary
{
public:
    FiniteSize();
    ~FiniteSize();

    // Boundary interface
public:
    int transformCoordinate(const int xi) const;
};

}
