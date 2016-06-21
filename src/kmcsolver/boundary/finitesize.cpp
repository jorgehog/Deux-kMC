#include "finitesize.h"

using namespace kMC;

FiniteSize::FiniteSize(const orientations orientation) :
    Boundary1D(orientation)
{

}

FiniteSize::~FiniteSize()
{

}

double FiniteSize::transformContinousCoordinate(const double xi) const
{
    return xi;
}

int FiniteSize::transformLatticeCoordinate(const int xi) const
{
    return xi;
}


void FiniteSize::closestImage(const double xi, const double xti, double &dxi) const
{
    noImage1D(xi, xti, dxi);
}
