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
    return std::fmod(xi + m_span, m_span);
}

bool Periodic::isBlocked(const double xi) const
{
    (void) xi;

    return false;
}

std::vector<double> Periodic::imagesOf(const double xi) const
{
    return {xi + (double)m_span, xi - (double)m_span};
}
