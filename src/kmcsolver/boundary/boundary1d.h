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
    virtual void closestImage(const double xi, const double xti, double &dxi) const = 0;

    void noImage1D(const double xi, const double xti, double &dxi) const;

    // Boundary interface
public:
    double transformCoordinate(const double xi, const double xj, const double xk) const;
    bool isBlocked(const double xi, const double xj, const double xk) const;
    void closestImage(const double xi, const double xj, const double xk,
                      const double xti, const double xtj, const double xtk,
                      double &dxi, double &dxj, double &dxk) const;
};

}
