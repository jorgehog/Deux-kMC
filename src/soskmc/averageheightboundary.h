#pragma once

#include "sosboundary.h"
#include "observers.h"

#include "subjects.h"

using kMC::Observer;
using kMC::Subjects;

class AverageHeightBoundary : public SOSBoundary, public Observer<Subjects>
{
public:
    AverageHeightBoundary(SOSSolver &solver, const uint averageHeightDepth,
                          const uint dim, const uint span, const uint yspan,
                          orientations orientation,
                          const uint location);
    virtual ~AverageHeightBoundary();

    template<typename T>
    const T &transformFunction(const T &xi, const T &xj, const T& xk) const;

    template<typename T>
    bool blockedFunction(const T& xi, const T &xj, const T& xk) const;

    void affectSurfaceSites();
    void affectSolutionSites(const int z);

    const double &average() const
    {
        return m_average;
    }

    double calcAverage() const;

    uint cutoffArea() const
    {
        return m_averageHeightDepth*m_yspan;
    }

    bool isInsideCutoff(const uint x, const uint y) const;

    void getStartsAndEnds(uint &x0,
                          uint &y0,
                          uint &x1,
                          uint &y1);

    const uint &dim() const
    {
        return m_dim;
    }

    const uint &location() const
    {
        return m_location;
    }

private:

    uint m_x0;
    uint m_y0;
    uint m_x1;
    uint m_y1;

    const uint m_averageHeightDepth;

    const uint m_dim;
    const uint m_span;
    const uint m_yspan;

    const uint m_location;
    double m_prevAverage;

    double m_average;

    // Boundary interface
public:
    double transformContinousCoordinate(const double xi, const double xj, const double xk) const;
    int transformLatticeCoordinate(const int xi, const int xj, const int xk) const;
    bool isBlockedContinous(const double xi, const double xj, const double xk) const;
    bool isBlockedLattice(const int xi, const int xj, const int xk) const;
    void closestImage(const double xi, const double xj, const double xk,
                      const double xti, const double xtj, const double xtk,
                      double &dxi, double &dxj, double &dxk) const;

    // kMC::Observer interface
public:
    void notifyObserver(const Subjects &subject);
    void initializeObserver(const Subjects &subject);
};

template<typename T>
const T &AverageHeightBoundary::transformFunction(const T &xi, const T &xj, const T& xk) const
{
    (void) xj;
    (void) xk;

    return xi;
}

template<typename T>
bool AverageHeightBoundary::blockedFunction(const T &xi, const T &xj, const T& xk) const
{
    (void) xj;

    const int averageHeight = round(average());

    if (orientation() == orientations::FIRST)
    {
        return (xi < T(m_location)) && (xk <= averageHeight);
    }

    else
    {
        return (xi > T(m_location)) && (xk <= averageHeight);
    }
}


