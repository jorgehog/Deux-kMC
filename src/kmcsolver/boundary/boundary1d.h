#pragma once

#include "boundary.h"


namespace kMC {

class Boundary1D : public Boundary
{
public:
    Boundary1D(const Boundary::orientations orientation);
    virtual ~Boundary1D();

    virtual double transformContinousCoordinate(const double xi) const = 0;
    virtual int transformLatticeCoordinate(const int xi) const = 0;
    virtual bool isBlockedContinous(const double xi) const = 0;
    virtual bool isBlockedLattice(const int xi) const = 0;
    virtual void closestImage(const double xi, const double xti, double &dxi) const = 0;

    void noImage1D(const double xi, const double xti, double &dxi) const;

    // Boundary interface
public:
    double transformContinousCoordinate(const double xi, const double xj, const double xk) const;
    int transformLatticeCoordinate(const int xi, const int xj, const int xk) const;
    bool isBlockedContinous(const double xi, const double xj, const double xk) const;
    bool isBlockedLattice(const int xi, const int xj, const int xk) const;
    void closestImage(const double xi, const double xj, const double xk,
                      const double xti, const double xtj, const double xtk,
                      double &dxi, double &dxj, double &dxk) const;
};

}
