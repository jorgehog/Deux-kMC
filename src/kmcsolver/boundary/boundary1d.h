#pragma once

#include "boundary.h"


namespace kMC {

class Boundary1D : public Boundary
{
public:
    Boundary1D(const Boundary::orientations orientation);
    virtual ~Boundary1D();

    virtual double transformCoordinate(const double xi) const = 0;
    virtual bool isBlocked(const double xi) const = 0;
    virtual std::vector<double> imagesOf(const double xi) const = 0;

    // Boundary interface
public:
    double transformCoordinate(const double xi, const double xj, const double xk) const;
    bool isBlocked(const double xi, const double xj, const double xk) const;
    std::vector<double> imagesOf(const double xi, const double xj, const double xk) const;
};

}
