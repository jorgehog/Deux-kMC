#pragma once

#include <sys/types.h>
#include <vector>
#include <limits>

namespace kMC
{

class KMCSolver;

class Boundary
{
public:

    enum class orientations
    {
        FIRST,
        LAST
    };

    Boundary(const Boundary::orientations orientation);

    virtual ~Boundary() {}

    virtual double transformContinousCoordinate(const double xi, const double xj, const double xk) const = 0;
    virtual int transformLatticeCoordinate(const int xi, const int xj, const int xk) const = 0;

    virtual bool isBlockedContinous(const double xi, const double xj, const double xk) const = 0;
    virtual bool isBlockedLattice(const int xi, const int xj, const int xk) const = 0;

    virtual void closestImage(const double xi, const double xj, const double xk,
                              const double xti, const double xtj, const double xtk,
                              double &dxi, double &dxj, double &dxk) const = 0;

    void noImage(const double xi, const double xj, const double xk,
                 const double xti, const double xtj, const double xtk,
                 double &dxi, double &dxj, double &dxk) const;

    const Boundary::orientations &orientation() const
    {
        return m_orientation;
    }

private:

    const Boundary::orientations m_orientation;

};

}
