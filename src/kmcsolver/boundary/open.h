#pragma once

#include "finitesize.h"

namespace kMC
{

class Open : public FiniteSize
{
public:
    Open(const Boundary::orientations orientation);
    ~Open();

    // Boundary interface
public:
    bool isBlocked(const double xi) const;
};

}
