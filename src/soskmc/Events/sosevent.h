#pragma once

#include "../sossolver.h"

#include <utils.h>

using namespace ignis;

class SOSEvent : public LatticeEvent
{
public:
    SOSEvent(const SOSSolver &solver,
                      string type = "Event",
                      string unit = "",
                      bool hasOutput=false,
                      bool storeValue=false) :
        LatticeEvent(type, unit, hasOutput, storeValue),
        m_solver(solver)

    {
        setDependency(solver);
    }

    virtual ~SOSEvent() {}

    const SOSSolver &solver() const
    {
        return m_solver;
    }

private:

    const SOSSolver &m_solver;

};


