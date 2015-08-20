#include "finitesize.h"

using namespace kMC;

FiniteSize::FiniteSize(const orientations orientation) :
    Boundary(orientation)
{

}

FiniteSize::~FiniteSize()
{

}

int FiniteSize::transformCoordinate(const int xi) const
{
    return xi;
}
