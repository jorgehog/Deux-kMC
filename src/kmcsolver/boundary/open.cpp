#include "open.h"

using namespace kMC;

Open::Open(const orientations orientation) :
    FiniteSize(orientation)
{

}

Open::~Open()
{

}

bool Open::isBlockedContinous(const double xi) const
{
    return blockedFunction<double>(xi);
}

bool Open::isBlockedLattice(const int xi) const
{
    return blockedFunction<int>(xi);
}
