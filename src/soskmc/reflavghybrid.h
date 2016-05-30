#pragma once


#include "averageheightboundary.h"

class ReflAvgHybrid : public AverageHeightBoundary
{
public:

    ReflAvgHybrid(SOSSolver &solver, const uint averageHeightDepth) :
        AverageHeightBoundary(solver, averageHeightDepth, 0,
                              solver.length(), solver.width(),
                              Boundary::orientations::LAST, solver.length() - 1)
    {

    }

    // Boundary interface
public:
    double transformContinousCoordinate(const double xi, const double xj, const double xk) const
    {
        (void) xj;
        (void) xk;

        if (xi >= solver().length() - 0.5)
        {

            return 2*solver().length() - xi - 1;
        }

        else
        {
            return xi;
        }
    }

    bool isBlockedContinous(const double xi, const double xj, const double xk) const
    {
        (void) xi;
        (void) xj;
        (void) xk;

        return false;
    }

    void closestImage(const double xi, const double xj, const double xk, const double xti, const double xtj, const double xtk, double &dxi, double &dxj, double &dxk) const
    {
        noImage(xi, xj, xk, xti, xtj, xtk, dxi, dxj, dxk);
    }
};
