#pragma once

#include "sosboundary.h"


class AverageHeightLineBoundary : public SOSBoundary, public Observer<Subjects>
{
public:
    AverageHeightLineBoundary(SOSSolver &solver, const Boundary::orientations orientation, const uint dim, const uint depth);

    const uint dim() const
    {
        return m_dim;
    }

private:

    const uint m_dim;
    const uint m_depth;

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

        bool isBelow = (xi <= m_location - 0.5) && (orientation() == orientations::FIRST);
        bool isAbove = (xi >= m_location - 0.5) && (orientation() == orientations::LAST);

        if (isBelow)
        {
            return 2*m_location - xi - 1;
        }

        else if (isAbove)
        {
            return 2*m_location - xi + 1;
        }

        else
        {
            return xi;
        }

        if (dim() == 0)
        {

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
};


