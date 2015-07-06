#include "confiningsurface.h"

ConfiningSurface::ConfiningSurface(SolidOnSolidSolver &solver,
                                   string type,
                                   string unit,
                                   bool hasOutput,
                                   bool storeValue) :
    SolidOnSolidEvent(solver, type, unit, hasOutput, storeValue),
    m_height(0)
{
    solver.setConfiningSurfaceEvent(*this);
}

ConfiningSurface::~ConfiningSurface()
{

}


