#pragma once


#include "boundary1d.h"

namespace kMC
{

class Periodic : public Boundary1D
{
public:

    Periodic(const uint span, const orientations orientation);

    ~Periodic();

    double transformCoordinate(const double xi) const;

    bool isBlocked(const double xi) const;

    void closestImage(const double xi, const double xti, double &dxi) const;

private:

    const uint m_span;

};

}

