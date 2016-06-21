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
    bool isBlockedContinous(const double xi) const
    {
        (void) xi;

        return false;
    }

    bool isBlockedLattice(const int xi) const
    {
        (void) xi;

        return false;
    }
};

}
