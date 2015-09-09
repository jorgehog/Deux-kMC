#pragma once


#include "boundary.h"

namespace kMC
{

class Periodic : public Boundary
{
public:

    Periodic(const uint span, const orientations orientation);

    ~Periodic();

    double transformCoordinate(const double xi) const;

    bool isBlocked(const double xi) const;

    std::vector<double> imagesOf(const double xi) const;

private:

    const uint m_span;

};

}

