#pragma once

#include "boundary.h"

namespace kMC
{


class ConstantHeight : public Boundary
{
public:
    ConstantHeight(const double height, const double location, const Boundary::orientations orientation);
    virtual ~ConstantHeight();

private:

    const double m_height;
    const double m_location;

    // Boundary interface
public:
    double transformCoordinate(const double xi, const double xj, const double xk) const;
    bool isBlocked(const double xi, const double xj, const double xk) const;
    void closestImage(const double xi, const double xj, const double xk,
                      const double xti, const double xtj, const double xtk,
                      double &dxi, double &dxj, double &dxk) const;
};

}

