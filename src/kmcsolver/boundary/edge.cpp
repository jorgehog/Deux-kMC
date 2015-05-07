#include "edge.h"

using namespace kMC;


Edge::Edge(const uint span) :
    Boundary(),
    m_span(span)
{

}

int Edge::transformCoordinate(const int xi) const
{
    return xi;
}

bool Edge::isBlocked(const int xi) const
{
    return xi < 0 || xi >= m_span;
}
