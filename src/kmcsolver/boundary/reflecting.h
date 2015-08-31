#pragma once

#include "boundary.h"

namespace kMC
{

class Reflecting : public Boundary
{
public:
    Reflecting(const int location, const Boundary::orientations orientation);
    ~Reflecting();

private:
    const int m_location;

    // Boundary interface
public:
    int transformCoordinate(const int xi) const;
    bool isBlocked(const int xi) const;
};

}

