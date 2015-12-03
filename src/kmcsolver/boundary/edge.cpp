#include "edge.h"

using namespace kMC;


Edge::Edge(const double location,
           const orientations orientation) :
    FiniteSize(orientation),
    m_location(location)
{

}



bool Edge::isBlockedContinous(const double xi) const
{
    return blockedFunction<double>(xi);
}

bool Edge::isBlockedLattice(const int xi) const
{
    return blockedFunction<int>(xi);
}
