#include "boundary1d.h"

using namespace kMC;


Boundary1D::Boundary1D(const Boundary::orientations orientation) :
    Boundary(orientation)
{

}

Boundary1D::~Boundary1D()
{

}


double kMC::Boundary1D::transformCoordinate(const double xi, const double xj, const double xk) const
{
    (void) xj;
    (void) xk;

    return transformCoordinate(xi);
}

bool kMC::Boundary1D::isBlocked(const double xi, const double xj, const double xk) const
{
    (void) xj;
    (void) xk;

    return isBlocked(xi);
}

std::vector<double> kMC::Boundary1D::imagesOf(const double xi, const double xj, const double xk) const
{
    (void) xj;
    (void) xk;

    return imagesOf(xi);
}
