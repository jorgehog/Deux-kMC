#pragma once

#include "../kmcsolver/reaction.h"

#include <vector>

using namespace kMC;

class SOSSolver;

class SOSReaction : public Reaction
{
public:
    SOSReaction(const uint x, const uint y, SOSSolver &system) :
        Reaction(),
        m_x(x),
        m_y(y),
        m_solver(system)
    {

    }

    SOSSolver &solver() const
    {
        return m_solver;
    }

    const uint &x() const
    {
        return m_x;
    }

    const uint &y() const
    {
        return m_y;
    }

    uint nNeighbors() const;

private:
    const uint m_x;
    const uint m_y;

    SOSSolver &m_solver;

};

