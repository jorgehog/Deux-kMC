#pragma once

#include <sys/types.h>

namespace kMC
{

class KMCSolver;

class Boundary
{
public:

    enum class orientations
    {
        FIRST,
        LAST
    };

    Boundary(const Boundary::orientations orientation);

    virtual ~Boundary() {}

    virtual int transformCoordinate(const int xi) const = 0;

    virtual bool isBlocked(const int xi) const = 0;


    const Boundary::orientations m_orientation;

};

}
