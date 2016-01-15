#include "confiningsurface.h"

#include "sossolver.h"

using namespace kMC;


ConfiningSurface::ConfiningSurface(SOSSolver &solver,
                                   string type,
                                   string unit,
                                   bool hasOutput,
                                   bool storeValue) :
    ignis::LatticeEvent(type, unit, hasOutput, storeValue),
    Observer(),
    Subject(),
    m_ccc(new CurrentConfinementChange()),
    m_height(0),
    m_solver(solver)
{
    solver.setConfiningSurfaceEvent(*this);
}

ConfiningSurface::~ConfiningSurface()
{
    delete m_ccc;
}

void ConfiningSurface::setHeight(const double height)
{
    m_ccc->prevHeight = height;

    m_height = height;

    notifyObservers(Subjects::CONFININGSURFACE);

    BADAss(solver().heights().max(), <=, height - 1);
}


