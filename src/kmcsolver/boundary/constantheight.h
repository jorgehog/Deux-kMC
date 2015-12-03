#pragma once

#include "boundary.h"


namespace kMC
{


class ConstantHeight : public Boundary
{
public:
    ConstantHeight(const double height, const double location, const Boundary::orientations orientation);
    virtual ~ConstantHeight();

    template<typename T>
    const T &transformFunction(const T &xi, const T &xj, const T &xk) const;

    template<typename T>
    bool blockedFunction(const T &xi, const T &xj, const T &xk) const;

private:

    const double m_height;
    const double m_location;

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

template<typename T>
bool ConstantHeight::blockedFunction(const T &xi, const T &xj, const T &xk) const
{
    (void) xj;

    if (orientation() == orientations::FIRST)
    {
        return (xi < m_location) && (xk <= m_height);
    }

    else
    {
        return (xi > m_location) && (xk <= m_height);
    }
}

template<typename T>
const T &ConstantHeight::transformFunction(const T &xi, const T &xj, const T &xk) const
{
    (void) xi;
    (void) xj;

    return xk;
}

}

