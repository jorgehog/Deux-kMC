#include "dumpsystem.h"

#include "diffusion/diffusion.h"

#include "sossolver.h"


using namespace kMC;


DumpSystem::DumpSystem(const SOSSolver &solver, const uint dumpInterval, const string path, const string ext) :
    SOSEvent(solver, "DumpSystem", "dumps", true),
    m_dumpInterval(dumpInterval),
    m_path(path),
    m_ext(ext)
{

}

DumpSystem::~DumpSystem()
{

}

uint DumpSystem::dumpInterval() const
{
    return m_dumpInterval;
}

void DumpSystem::dump(const uint n)
{
    solver().diffusionEvent().dump(n, m_path, m_ext);
    m_nDumps++;
}

void DumpSystem::initialize()
{
    m_nDumps = 0;
    dump(0);
}

void DumpSystem::execute()
{
    setValue(m_nDumps);
}

void DumpSystem::reset()
{
    if ((cycle()+1) % dumpInterval() == 0)
    {
        dump((cycle()+1)/dumpInterval());
    }
}
