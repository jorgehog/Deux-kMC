#pragma once

#include "finitesize.h"

namespace kMC
{

class Open : public FiniteSize
{
public:
    Open(const Boundary::orientations orientation);
    ~Open();

    template<typename T>
    bool blockedFunction(const T& xi) const;

    // Boundary interface
public:
    bool isBlockedContinous(const double xi) const;
    bool isBlockedLattice(const int xi) const;

};

template<typename T>
bool Open::blockedFunction(const T &xi) const
{
    (void) xi;

    return false;
}


}
