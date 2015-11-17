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
                          const double location);
    virtual ~AverageHeightBoundary();

    void affectSurfaceSites();
    void affectSolutionSites(const int z);
    void updateSites(const int height, const int prevHeight);

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

private:

    uint m_x0;
    uint m_y0;
    uint m_x1;
    uint m_y1;

    SOSSolver &m_mutexSolver;

    const uint m_averageHeightDepth;

    const uint m_dim;
    const uint m_span;
    const uint m_yspan;

    const double m_location;
    double m_prevAverage;

    double m_average;

    // Boundary interface
public:
    double transformCoordinate(const double xi, const double xj, const double xk) const;
    bool isBlocked(const double xi, const double xj, const double xk) const;
    void closestImage(const double xi, const double xj, const double xk,
                      const double xti, const double xtj, const double xtk,
                      double &dxi, double &dxj, double &dxk) const;

    // kMC::Observer interface
public:
    void notifyObserver(const Subjects &subject);
    void initializeObserver(const Subjects &subject);
};

