#pragma once

#include "sosevent.h"

namespace kMC
{

class DumpSystem : public SOSEvent
{
public:
    DumpSystem(const SOSSolver &solver, const uint dumpInterval, const string path = "/tmp");
    ~DumpSystem();

    uint dumpInterval() const;

private:

    const uint m_dumpInterval;
    uint m_nDumps;

    const string m_path;

    // Event interface
public:
    void initialize();
    void execute();
};

}

