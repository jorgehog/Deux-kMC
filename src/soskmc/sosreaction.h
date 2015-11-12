#pragma once

#include "../kmcsolver/reaction.h"

#include <vector>

using namespace kMC;

class SOSSolver;

class SOSReaction : public Reaction
{
public:
    SOSReaction(const uint x, const uint y, SOSSolver &solver) :
        Reaction(),
        m_x(x),
        m_y(y),
        m_solver(solver)
    {

    }

    virtual ~SOSReaction();

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

    void setX(const uint x)
    {
        m_x = x;
    }

    void setY(const uint y)
    {
        m_y = y;
    }

private:
    uint m_x;
    uint m_y;

    SOSSolver &m_solver;

};

