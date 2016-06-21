#pragma once

#include "boundary1d.h"

namespace kMC
{

class Reflecting : public Boundary1D
{
public:
    Reflecting(const int location, const Boundary::orientations orientation);
    ~Reflecting();

    const int &location() const
    {
        return m_location;
    }

private:
    const int m_location;

    // Boundary interface
public:
    double transformContinousCoordinate(const double xi) const;
    int transformLatticeCoordinate(const int xi) const;
    bool isBlockedContinous(const double xi) const;
    bool isBlockedLattice(const int xi) const;
    void closestImage(const double xi, const double xti, double &dxi) const;
};
}

