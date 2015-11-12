#include "confiningsurface.h"

#include "sossolver.h"


ConfiningSurface::ConfiningSurface(SOSSolver &solver,
                                   string type,
                                   string unit,
                                   bool hasOutput,
                                   bool storeValue) :
    ignis::LatticeEvent(type, unit, hasOutput, storeValue),
    m_height(0),
    m_solver(solver)
{
    solver.setConfiningSurfaceEvent(*this);
}

ConfiningSurface::~ConfiningSurface()
{

}


