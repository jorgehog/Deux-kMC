#include "sosreaction.h"

#include "sossolver.h"

#include "../kmcsolver/boundary/boundary.h"


SOSReaction::~SOSReaction()
{

}

uint SOSReaction::nNeighbors() const
{
    return m_solver.nNeighbors(m_x, m_y);
}
