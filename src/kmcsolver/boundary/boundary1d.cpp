#include "boundary1d.h"

using namespace kMC;


Boundary1D::Boundary1D(const Boundary::orientations orientation) :
    Boundary(orientation)
{

}

Boundary1D::~Boundary1D()
{

}

void Boundary1D::noImage1D(const double xi, const double xti, double &dxi) const
{
    dxi = xti - xi;
}


double kMC::Boundary1D::transformContinousCoordinate(const double xi, const double xj, const double xk) const
{
    (void) xj;
    (void) xk;

    return transformContinousCoordinate(xi);
}

int Boundary1D::transformLatticeCoordinate(const int xi, const int xj, const int xk) const
{
    (void) xj;
    (void) xk;

    return transformLatticeCoordinate(xi);
}

bool Boundary1D::isBlockedContinous(const double xi, const double xj, const double xk) const
{
    (void) xj;
    (void) xk;

    return isBlockedContinous(xi);
}

bool Boundary1D::isBlockedLattice(const int xi, const int xj, const int xk) const
{
    (void) xj;
    (void) xk;

    return isBlockedLattice(xi);
}

void Boundary1D::closestImage(const double xi, const double xj, const double xk,
                              const double xti, const double xtj, const double xtk,
                              double &dxi, double &dxj, double &dxk) const
{
    (void) xj;
    (void) xk;
    (void) dxj;
    (void) dxk;
    (void) xtj;
    (void) xtk;

    dxj = std::numeric_limits<double>::max();
    dxk = std::numeric_limits<double>::max();

    closestImage(xi, xti, dxi);
}
