#pragma once

#include "boundary1d.h"

namespace kMC
{

class Reflecting : public Boundary1D
{
public:
    Reflecting(const int location, const Boundary::orientations orientation);
    ~Reflecting();

    template<typename T>
    const T transformFunction(const T &xi) const;

    template<typename T>
    bool blockedFunction(const T& xi) const;


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

template<typename T>
const T Reflecting::transformFunction(const T &xi) const
{
    bool isBelow = (xi < m_location) && (orientation() == orientations::FIRST);
    bool isAbove = (xi > m_location) && (orientation() == orientations::LAST);

    if (isBelow)
    {
        return 2*m_location - xi - 1;
    }

    else if (isAbove)
    {
        return 2*m_location - xi + 1;
    }

    else
    {
        return xi;
    }
}

template<typename T>
bool Reflecting::blockedFunction(const T &xi) const
{
    (void) xi;

    return false;
}

}

