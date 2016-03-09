#pragma once

#include "sosboundary.h"


class AverageHeightLineBoundary : public SOSBoundary
{
public:
    AverageHeightLineBoundary(SOSSolver &solver, const Boundary::orientations orientation, const uint dim, const uint depth);

private:

    const uint m_dim;
    const uint m_depth;

    uint m_x0;
    uint m_x1;

    // Boundary interface
public:
    double transformContinousCoordinate(const double xi, const double xj, const double xk) const;
    int transformLatticeCoordinate(const int xi, const int xj, const int xk) const;
    bool isBlockedContinous(const double xi, const double xj, const double xk) const;
    bool isBlockedLattice(const int xi, const int xj, const int xk) const;
    void closestImage(const double xi, const double xj, const double xk, const double xti, const double xtj, const double xtk, double &dxi, double &dxj, double &dxk) const;
};
