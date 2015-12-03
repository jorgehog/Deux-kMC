#pragma once

#include "boundary1d.h"

namespace kMC
{

class FiniteSize : public Boundary1D
{
public:
    FiniteSize(const Boundary::orientations orientation);
    ~FiniteSize();

    template<typename T>
    const T& transformFunction(const T& xi) const;

    // Boundary interface
public:
    double transformContinousCoordinate(const double xi) const;
    int transformLatticeCoordinate(const int xi) const;
    void closestImage(const double xi, const double xti, double &dxi) const;
};


template<typename T>
const T &FiniteSize::transformFunction(const T &xi) const
{
    return xi;
}


}
