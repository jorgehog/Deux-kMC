#pragma once


#include "boundary1d.h"

namespace kMC
{

class Periodic : public Boundary1D
{
public:

    Periodic(const int span, const orientations orientation);

    ~Periodic();

    template<typename T>
    bool blockedFunction(const T& xi) const;

    double transformContinousCoordinate(const double xi) const;

    int transformLatticeCoordinate(const int xi) const;

    bool isBlockedContinous(const double xi) const;

    bool isBlockedLattice(const int xi) const;

    void closestImage(const double xi, const double xti, double &dxi) const;

private:

    const int m_span;
    const double m_spanMinusHalf;

};


template<typename T>
bool Periodic::blockedFunction(const T &xi) const
{
    (void) xi;

    return false;
}

}

