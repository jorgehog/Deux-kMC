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



double Reflecting::transformCoordinate(const double xi) const
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

bool Reflecting::isBlocked(const double xi) const
{
    (void) xi;

    return false;
}

std::vector<double> Reflecting::imagesOf(const double xi) const
{
    if (m_orientation == orientations::FIRST)
    {
        return {2*m_location - xi - 1};
    }

    else
    {
        return {2*m_location - xi + 1};
    }

}
