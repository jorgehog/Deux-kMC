#include "periodic.h"

#include <cmath>

using namespace kMC;

Periodic::Periodic(const int span, const Boundary::orientations orientation) :
    Boundary1D(orientation),
    m_span(span),
    m_spanMinusHalf(span-0.5)
{

}

Periodic::~Periodic()
{

}

double Periodic::transformContinousCoordinate(const double xi) const
{
    if (xi > m_spanMinusHalf)
    {
        return xi - m_span;
    }

    else if (xi < -0.5)
    {
        return xi + m_span;
    }

    else
    {
        return xi;
    }
}

int Periodic::transformLatticeCoordinate(const int xi) const
{
    if (xi >= m_span)
    {
        return xi - m_span;
    }

    else if (xi < 0)
    {
        return xi + m_span;
    }

    else
    {
        return xi;
    }

//    return (xi + m_span) % m_span;
}

bool Periodic::isBlockedContinous(const double xi) const
{
    return blockedFunction(xi);
}

bool Periodic::isBlockedLattice(const int xi) const
{
    return blockedFunction(xi);
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
