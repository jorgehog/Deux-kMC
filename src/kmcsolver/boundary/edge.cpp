#include "edge.h"

using namespace kMC;


Edge::Edge(const double location,
           const orientations orientation) :
    FiniteSize(orientation),
    m_location(location)
{

}

bool Edge::isBlocked(const double xi) const
{
    if (orientation() == orientations::FIRST)
    {
        return xi < m_location;
    }

    else
    {
        return xi > m_location;
    }
}
