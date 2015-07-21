#pragma once

#include "../sossolver.h"

#include <utils.h>

using namespace ignis;

class SolidOnSolidEvent : public LatticeEvent
{
public:
    SolidOnSolidEvent(const SOSSolver &solver,
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

    const SOSSolver &solver() const
    {
        return m_solver;
    }

private:

    const SOSSolver &m_solver;

};

