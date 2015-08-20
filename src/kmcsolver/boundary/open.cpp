#include "open.h"

using namespace kMC;

Open::Open(const orientations orientation) :
    FiniteSize(orientation)
{

}

Open::~Open()
{

}



bool Open::isBlocked(const int xi) const
{
    (void) xi;

    return false;
}
