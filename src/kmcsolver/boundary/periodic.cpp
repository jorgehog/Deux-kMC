#include "periodic.h"

#include <cmath>

using namespace kMC;

Periodic::Periodic(const uint span, const Boundary::orientations orientation) :
    Boundary1D(orientation),
    m_span(span)
{

}

Periodic::~Periodic()
{

}

double Periodic::transformCoordinate(const double xi) const
{
    return std::fmod(xi + m_span + 0.5, m_span) - 0.5;
}

bool Periodic::isBlocked(const double xi) const
{
    (void) xi;

    return false;
}

void Periodic::closestImage(const double xi, const double xti, double &dxi) const
{
    double dxti = xti - xi;

    if (fabs(dxti) < m_span/2)
    {
        dxi = dxti;
    }

    else
    {
        if (xi < m_span/2)
        {
            dxi = xti - (double)m_span - xi;
        }

        else
        {
            dxi = xti + (double)m_span - xi;
        }
    }
}
