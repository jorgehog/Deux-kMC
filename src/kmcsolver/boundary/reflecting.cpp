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
    if (m_location == 0)
    {
        if (xi < -0.5)
        {
            return -(xi + 1);
        }
    }

    else
    {
        if (xi > m_location + 0.5)
        {
            return 2*m_location - (xi - 1);
        }
    }

    return xi;
}

int Reflecting::transformLatticeCoordinate(const int xi) const
{
    if (m_location == 0)
    {
        if (xi < 0)
        {
            return -(xi + 1);
        }
    }

    else
    {
        if (xi > m_location)
        {
            return 2*m_location - (xi - 1);
        }
    }

    return xi;
}

bool Reflecting::isBlockedContinous(const double xi) const
{
    (void) xi;
    return false;
}

bool Reflecting::isBlockedLattice(const int xi) const
{
    (void) xi;
    return false;
}



void Reflecting::closestImage(const double xi, const double xti, double &dxi) const
{
    noImage1D(xi,xti, dxi);
}
