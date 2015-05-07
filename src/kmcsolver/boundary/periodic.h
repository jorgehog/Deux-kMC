#pragma once


#include "boundary.h"

namespace kMC
{

class Periodic : public Boundary
{
public:

    Periodic(const uint span);

    ~Periodic() {}

    int transformCoordinate(const int xi) const;

    bool isBlocked(const int xi) const;

private:

    const uint m_span;


};

}

