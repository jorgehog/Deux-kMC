#include "reflecting.h"

using namespace kMC;

Reflecting::Reflecting(const int location, const Boundary::orientations orientation) :
    Boundary(orientation),
    m_location(location)
{

}

Reflecting::~Reflecting()
{

}



int kMC::Reflecting::transformCoordinate(const int xi) const
{
    bool isBelow = (xi < m_location) && (m_orientation == orientations::FIRST);
    bool isAbove = (xi > m_location) && (m_orientation == orientations::LAST);

    if (isBelow)
    {
        return 2*m_location - xi - 1;
    }

    else if (isAbove)
    {
        return 2*m_location - xi + 1;
    }

    else
    {
        return xi;
    }
}

bool kMC::Reflecting::isBlocked(const int xi) const
{
    (void) xi;

    return false;
}
