#pragma once

#include "sosboundary.h"
#include "averageheightboundary.h"


class AverageHeightLineBoundary : public SOSBoundary, public Observer<Subjects>
{
public:
    AverageHeightLineBoundary(SOSSolver &solver, const Boundary::orientations orientation, const uint dim, const uint depth);

    const uint &dim() const
    {
        return m_dim;
    }

    const uint &location() const
    {
        return m_location;
    }

private:

    const uint m_dim;
    const uint m_depth;

    uint m_location;

    uint m_x0;
    uint m_x1;

    // Boundary interface
public:
    int transformLatticeCoordinate(const int xi, const int xj, const int xk) const;
    bool isBlockedLattice(const int xi, const int xj, const int xk) const;
    void closestImage(const double xi, const double xj, const double xk, const double xti, const double xtj, const double xtk, double &dxi, double &dxj, double &dxk) const;

    // Observer interface
public:
    void initializeObserver(const Subjects &subject);
    void notifyObserver(const Subjects &subject);
};

class AverageHeightLineBoundaryOpen : public AverageHeightLineBoundary
{
public:
    using AverageHeightLineBoundary::AverageHeightLineBoundary;

    // Boundary interface
    double transformContinousCoordinate(const double xi, const double xj, const double xk) const
    {
        (void) xj;
        (void) xk;
        return xi;
    }


    bool isBlockedContinous(const double xi, const double xj, const double xk) const
    {
        (void) xi;
        (void) xj;
        (void) xk;

        return false;
    }
};

class AverageHeightLineBoundaryRefl : public AverageHeightLineBoundary
{
public:
    using AverageHeightLineBoundary::AverageHeightLineBoundary;

    // Boundary interface
    double transformContinousCoordinate(const double xi, const double xj, const double xk) const
    {
        (void) xj;
        (void) xk;

        if (location() == 0)
        {
            if (xi < -0.5)
            {
                return -(xi + 1);
            }
        }

        else
        {
            if (xi > location() - 0.5)
            {
                return 2*location() - (xi + 1);
            }
        }

        return xi;
    }

    bool isBlockedContinous(const double xi, const double xj, const double xk) const
    {
        (void) xi;
        (void) xj;
        (void) xk;

        return false;
    }
};


class ReflConstantHybrid : public SOSBoundary
{
public:

    ReflConstantHybrid(SOSSolver &solver, const int h, const Boundary::orientations orientation, const uint dim);

    const uint &location() const
    {
        return m_location;
    }

private:
    const int m_h;
    uint m_location;

    // Boundary interface
public:
    double transformContinousCoordinate(const double xi, const double xj, const double xk) const
    {
        (void) xj;
        (void) xk;

        if (location() == 0)
        {
            if (xi < -0.5)
            {
                return -(xi + 1);
            }
        }

        else
        {
            if (xi > location() - 0.5)
            {
                return 2*location() - (xi + 1);
            }
        }

        return xi;
    }

    int transformLatticeCoordinate(const int xi, const int xj, const int xk) const
    {
        (void) xj;
        (void) xk;

        return xi;
    }

    bool isBlockedContinous(const double xi, const double xj, const double xk) const
    {
        (void) xi;
        (void) xj;
        (void) xk;

        return false;
    }

    bool isBlockedLattice(const int xi, const int xj, const int xk) const
    {
        (void) xj;


        return (xi >= (int)location()) && (xk <= m_h);
    }

    void closestImage(const double xi, const double xj, const double xk, const double xti, const double xtj, const double xtk, double &dxi, double &dxj, double &dxk) const
    {
        noImage(xi, xj, xk, xti, xtj, xtk, dxi, dxj, dxk);
    }
};


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

    bool isBlockedLattice(const int xi, const int xj, const int xk) const
    {
        (void) xj;
        (void) xk;

        return xi >= (int)solver().length();
    }

    void closestImage(const double xi, const double xj, const double xk, const double xti, const double xtj, const double xtk, double &dxi, double &dxj, double &dxk) const
    {
        noImage(xi, xj, xk, xti, xtj, xtk, dxi, dxj, dxk);
    }
};
