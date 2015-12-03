#pragma once

#include "finitesize.h"

namespace kMC
{

class Edge : public FiniteSize
{
public:

    Edge(const double location,
         const orientations orientation);

    template<typename T>
    bool blockedFunction(const T& xi) const;

private:

    const double m_location;


    // Boundary1D interface
public:
    bool isBlockedContinous(const double xi) const;
    bool isBlockedLattice(const int xi) const;
};

template<typename T>
bool Edge::blockedFunction(const T &xi) const
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

}
