#pragma once

#include "Events/sosevent.h"

class SOSSolver;

class ParticleNumberConservator : public LatticeEvent
{
public:
    ParticleNumberConservator(SOSSolver &solver);

private:

    SOSSolver &m_solver;
    uint m_targetN;

    SOSSolver &solver() const
    {
        return m_solver;
    }

    // Event interface
public:
    void execute();
    void initialize();
    void reset();
};
