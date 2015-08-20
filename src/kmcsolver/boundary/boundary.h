#pragma once

#include <sys/types.h>

namespace kMC
{

class KMCSolver;

class Boundary
{
public:

    Boundary();

    virtual ~Boundary() {}

    virtual int transformCoordinate(const int xi) const = 0;

    virtual bool isBlocked(const int xi) const = 0;

    enum class orientations
    {
        FIRST,
        LAST
    };

};

}
