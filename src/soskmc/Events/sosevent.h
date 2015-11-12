#pragma once

#include <utils.h>

using namespace ignis;

class SOSSolver;

class SOSEvent : public LatticeEvent
{
public:
    SOSEvent(const SOSSolver &solver,
             string type = "Event",
             string unit = "",
             bool hasOutput=false,
             bool storeValue=false);

    virtual ~SOSEvent() {}

    const SOSSolver &solver() const
    {
        return m_solver;
    }

private:

    const SOSSolver &m_solver;

};

class MutexSOSEvent : public SOSEvent
{
public:

    MutexSOSEvent(SOSSolver &solver,
                  string type = "Event",
                  string unit = "",
                  bool hasOutput=false,
                  bool storeValue=false);

    SOSSolver &solver() const
    {
        return m_mutexSolver;
    }

private:

    SOSSolver &m_mutexSolver;
};

