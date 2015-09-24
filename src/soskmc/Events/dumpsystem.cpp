#include "dumpsystem.h"

#include "diffusion/diffusion.h"

using namespace kMC;


DumpSystem::DumpSystem(const SOSSolver &solver, const uint dumpInterval, const string path) :
    SOSEvent(solver, "DumpSystem", "dumps", true),
    m_dumpInterval(dumpInterval),
    m_path(path)
{

}

DumpSystem::~DumpSystem()
{

}

uint DumpSystem::dumpInterval() const
{
    return m_dumpInterval;
}

void DumpSystem::initialize()
{
    m_nDumps = 0;
}

void kMC::DumpSystem::execute()
{
    if (cycle() % dumpInterval() == 0)
    {
        solver().diffusionEvent().dump(cycle()/dumpInterval(), m_path);
        m_nDumps++;
    }

    setValue(m_nDumps);
}
