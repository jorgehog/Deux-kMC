#pragma once

#include "finitesize.h"

namespace kMC
{

class Open : public FiniteSize
{
public:
    Open();
    ~Open();

    // Boundary interface
public:
    bool isBlocked(const int xi) const;
};

}
