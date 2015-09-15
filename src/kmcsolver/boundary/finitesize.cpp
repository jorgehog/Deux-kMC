#include "finitesize.h"

using namespace kMC;

FiniteSize::FiniteSize(const orientations orientation) :
    Boundary1D(orientation)
{

}

FiniteSize::~FiniteSize()
{

}

double FiniteSize::transformCoordinate(const double xi) const
{
    return xi;
}

std::vector<double> FiniteSize::imagesOf(const double xi) const
{
    (void) xi;

    return {};
}
