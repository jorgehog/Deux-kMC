#pragma once

#include "../solidonsolidsolver.h"

#include <utils.h>

using namespace ignis;

class SolidOnSolidEvent : public LatticeEvent
{
public:
    SolidOnSolidEvent(const SolidOnSolidSolver &solver,
                      string type = "Event",
                      string unit = "",
                      bool hasOutput=false,
                      bool storeValue=false) :
        LatticeEvent(type, unit, hasOutput, storeValue),
        m_solver(solver)

    {
        setDependency(solver);
    }

    virtual ~SolidOnSolidEvent() {}

    const SolidOnSolidSolver &solver() const
    {
        return m_solver;
    }

private:

    const SolidOnSolidSolver &m_solver;

};
