#pragma once

#include "../kmcsolver/boundary/reflecting.h"
#include "observers.h"
#include "subjects.h"
#include "localpotential.h"

#include <BADAss/badass.h>

#include <armadillo>
using arma::vec;
using arma::ivec;

namespace kMC
{

class BoundaryTrackingDevice : public Observer<Subjects>
{
public:
    BoundaryTrackingDevice(SOSSolver &solver,
                           const uint location,
                           const uint dim,
                           const uint depth);

    virtual double averageHeight(const uint boundarySite) const = 0;

    double averageHeight(const uint x, const uint y) const
    {
        if (dim() == 0)
        {
            return averageHeight(y);
        }

        else
        {
            return averageHeight(x);
        }
    }

    bool pointIsInsideArea(const int x, const int y) const;

    bool pointIsAtBoundary(const uint x, const uint y) const
    {
        if (dim() == 0)
        {
            return x == location();
        }

        else
        {
            return y == location();
        }
    }

    void onAverageChange(const bool affect);

    const uint &area() const
    {
        return m_area;
    }

    const uint &dim() const
    {
        return m_dim;
    }

    const uint &location() const
    {
        return m_location;
    }

    uint orientation() const
    {
        if (m_location == 0)
        {
            return 0;
        }

        else
        {
            return 1;
        }
    }

    const SOSSolver &solver() const
    {
        return m_solver;
    }

    const uint &x0() const
    {
        return m_x0;
    }

    const uint &x1() const
    {
        return m_x1;
    }

    const uint &depth() const
    {
        return m_depth;
    }

protected:
    SOSSolver &solver()
    {
        return m_solver;
    }

private:
    SOSSolver &m_solver;
    const uint m_location;
    const uint m_dim;
    const uint m_depth;
    const uint m_area;
    uint m_x0;
    uint m_x1;
};

class TrackLineAverage : public BoundaryTrackingDevice
{
public:
    TrackLineAverage(SOSSolver &solver,
                     const uint location,
                     const uint dim,
                     const uint depth);

    double bruteForceAverage(const uint boundarySite) const;

private:
    vec m_averages;

    // Observer interface
public:
    void initializeObserver(const Subjects &subject);

    void notifyObserver(const Subjects &subject);

    // TrackingDevice interface
public:
    double averageHeight(const uint boundarySite) const
    {
        BADAssClose(m_averages(boundarySite), bruteForceAverage(boundarySite), 1E-10);
        return m_averages(boundarySite);
    }
};

class TrackAreaAverage : public BoundaryTrackingDevice
{
public:
    TrackAreaAverage(SOSSolver &solver,
                     const uint location,
                     const uint dim,
                     const uint depth) :
        BoundaryTrackingDevice(solver, location, dim, depth)
    {

    }

    double getBruteForceAverage() const;

private:
    double m_average;

    // Observer interface
public:
    void initializeObserver(const Subjects &subject)
    {
        (void) subject;

        m_average = getBruteForceAverage();

        onAverageChange(false);

    }
    void notifyObserver(const Subjects &subject);

    // TrackingDevice interface
public:
    double averageHeight(const uint boundarySite) const
    {
        (void) boundarySite;

        BADAssClose(m_average, getBruteForceAverage(), 1E-10);
        return m_average;
    }
};


//!This class removes a neighbor bond across a blocked boundary.
class NoBoundaryNeighbors : public LocalPotential
{
public:
    NoBoundaryNeighbors(SOSSolver &solver,
                        const int surfaceLevel,
                        const uint location,
                        const uint dim) :
        LocalPotential(solver),
        m_surfaceLevel(surfaceLevel),
        m_location(location),
        m_dim(dim)
    {

    }

    bool pointIsOnBoundary(const uint x, const uint y) const
    {
        if (m_dim == 0)
        {
            return x == m_location;
        }

        else
        {
            return y == m_location;
        }
    }

private:

    const int m_surfaceLevel;
    const uint m_location;
    const uint m_dim;

    // LocalPotential interface
public:
    double potential(const uint x, const uint y) const;
};


//!See comments in potential function for explanation.
class PartialBoundaryNeighbors : public LocalPotential
{
public:
    PartialBoundaryNeighbors(SOSSolver &solver, const BoundaryTrackingDevice &tracker);

    uint selectBoundarySite(const uint x, const uint y) const;

private:
    const BoundaryTrackingDevice &m_tracker;
    ivec m_boundaryNeighbors;

    // LocalPotential interface
public:
    double potential(const uint x, const uint y) const;
};


class ReflectingSurfaceOpenSolution : public Reflecting
{
    using Reflecting::Reflecting;

    // Boundary interface
public:
    double transformContinousCoordinate(const double xi) const
    {
        return xi;
    }
    bool isBlockedContinous(const double xi) const
    {
        (void) xi;
        return false;
    }
};


class BlockByTracker : public Boundary
{
public:
    BlockByTracker(const BoundaryTrackingDevice &tracker,
                   const double treshold = 0.5) :
        Boundary(tracker.location() == 0 ? First : Last),
        m_tracker(tracker),
        m_treshold(treshold)
    {

    }

private:
    const BoundaryTrackingDevice &m_tracker;
    const double m_treshold;


    // Boundary interface
public:
    double transformContinousCoordinate(const double xi, const double xj, const double xk) const
    {
        (void) xj;
        (void) xk;

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

    bool isBlockedLattice(const int xi, const int xj, const int xk) const;

    void closestImage(const double xi, const double xj, const double xk, const double xti, const double xtj, const double xtk, double &dxi, double &dxj, double &dxk) const
    {
        return noImage(xi, xj, xk, xti, xtj, xtk, dxi, dxj, dxk);
    }
};

}
