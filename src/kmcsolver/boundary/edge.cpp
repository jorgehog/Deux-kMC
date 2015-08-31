#include "edge.h"

using namespace kMC;


Edge::Edge(const int location,
           const orientations orientation) :
    FiniteSize(orientation),
    m_location(location)
{

}

bool Edge::isBlocked(const int xi) const
{
    if (m_orientation == orientations::FIRST)
    {
        return xi < m_location;
    }

    else
    {
        return xi > m_location;
    }
}
