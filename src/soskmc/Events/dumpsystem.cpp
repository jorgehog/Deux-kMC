#include "dumpsystem.h"

#include "diffusion/diffusion.h"

using namespace kMC;


DumpSystem::DumpSystem(const SOSSolver &solver, const uint dumpInterval) :
    SOSEvent(solver, "DumpSystem"),
    m_dumpInterval(dumpInterval)
{

}

DumpSystem::~DumpSystem()
{

}

uint DumpSystem::dumpInterval() const
{
    return m_dumpInterval;
}

void kMC::DumpSystem::execute()
{
    if (cycle() % dumpInterval() == 0)
    {
        solver().diffusionEvent().dump(cycle()/dumpInterval());
    }
}
