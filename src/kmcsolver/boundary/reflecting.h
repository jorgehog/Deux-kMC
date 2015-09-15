#pragma once

#include "boundary1d.h"

namespace kMC
{

class Reflecting : public Boundary1D
{
public:
    Reflecting(const int location, const Boundary::orientations orientation);
    ~Reflecting();

private:
    const int m_location;

    // Boundary interface
public:
    double transformCoordinate(const double xi) const;
    bool isBlocked(const double xi) const;
    std::vector<double> imagesOf(const double xi) const;
};

}

