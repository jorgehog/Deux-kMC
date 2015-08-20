#include "periodic.h"

using namespace kMC;

Periodic::Periodic(const uint span, const Boundary::orientations orientation) :
    Boundary(orientation),
    m_span(span)
{

}

int Periodic::transformCoordinate(const int xi) const
{
    return (xi + m_span)%m_span;
}

bool Periodic::isBlocked(const int xi) const
{
    (void) xi;

    return false;
}
