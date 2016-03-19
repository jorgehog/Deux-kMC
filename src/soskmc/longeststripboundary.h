#pragma once

#include "sosboundary.h"
#include "observers.h"

#include "subjects.h"

using kMC::Observer;
using kMC::Subjects;


class LongestStripBoundary : public SOSBoundary, public Observer<Subjects>
{
public:
    LongestStripBoundary(SOSSolver &solver,
                         const uint dim,
                         orientations orientation);

    bool isConnectedToOriginCalc(const uint &y, const int &h) const;
    bool isConnectedToOrigin(const uint y, const int h) const;

private:
    double m_location;
    const uint m_dim;
    uint m_yspan;
    uint m_ymax;
    uint m_ymax2;
    int m_currentHeight;

    void scanForYMax();
    void scanForYMax2();

    // Boundary interface
public:
    double transformContinousCoordinate(const double xi, const double xj, const double xk) const;
    int transformLatticeCoordinate(const int xi, const int xj, const int xk) const;
    bool isBlockedContinous(const double xi, const double xj, const double xk) const;
    bool isBlockedLattice(const int xi, const int xj, const int xk) const;
    void closestImage(const double xi, const double xj, const double xk, const double xti, const double xtj, const double xtk, double &dxi, double &dxj, double &dxk) const;

    // Observer interface
public:
    void initializeObserver(const Subjects &subject);
    void notifyObserver(const Subjects &subject);
};

