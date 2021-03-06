#pragma once

#include "sosevent.h"

namespace kMC
{

class DumpSystem : public SOSEvent
{
public:
    DumpSystem(const SOSSolver &solver, const uint dumpInterval, const string path = "/tmp", const string ext = "");
    ~DumpSystem();

    uint dumpInterval() const;

private:

    const uint m_dumpInterval;
    uint m_nDumps;

    const string m_path;
    const string m_ext;

    void dump(const uint n);

    // Event interface
public:
    void initialize();
    void execute();
    void reset();
};

}

