#pragma once

#include <sys/types.h>
#include <vector>

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

    virtual double transformCoordinate(const double xi, const double xj, const double xk) const = 0;

    virtual bool isBlocked(const double xi, const double xj, const double xk) const = 0;

    virtual std::vector<double> imagesOf(const double xi, const double xj, const double xk) const = 0;


    const Boundary::orientations m_orientation;

};

}
