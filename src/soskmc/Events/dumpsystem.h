#pragma once

#include "sosevent.h"

namespace kMC
{

class DumpSystem : public SOSEvent
{
public:
    DumpSystem(const SOSSolver &solver, const uint dumpInterval);
    ~DumpSystem();

    uint dumpInterval() const;

private:

    const uint m_dumpInterval;


    // Event interface
public:
    void execute();
};

}

