#include "reflecting.h"

using namespace kMC;

Reflecting::Reflecting(const int location, const Boundary::orientations orientation) :
    Boundary1D(orientation),
    m_location(location)
{

}

Reflecting::~Reflecting()
{

}

double Reflecting::transformContinousCoordinate(const double xi) const
{
    return transformFunction(xi);
}

int Reflecting::transformLatticeCoordinate(const int xi) const
{
    return transformFunction(xi);
}

bool Reflecting::isBlockedContinous(const double xi) const
{
    return blockedFunction(xi);
}

bool Reflecting::isBlockedLattice(const int xi) const
{
    return blockedFunction(xi);
}



void Reflecting::closestImage(const double xi, const double xti, double &dxi) const
{
    noImage1D(xi,xti, dxi);
}
