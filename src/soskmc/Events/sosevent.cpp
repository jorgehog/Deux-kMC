#include "sosevent.h"

#include "../sossolver.h"
#include "../Events/diffusion/diffusion.h"
#include "../Events/confiningsurface/confiningsurface.h"

SOSEvent::SOSEvent(const SOSSolver &solver,
                   std::string type,
                   std::string unit,
                   bool hasOutput,
                   bool storeValue) :
    LatticeEvent(type, unit, hasOutput, storeValue),
    m_solver(solver)

{
    BADAss(&solver.confiningSurfaceEvent(), !=, nullptr);
    BADAss(&solver.diffusionEvent(), !=, nullptr);

    setDependency(solver);
    setDependency(solver.confiningSurfaceEvent());
    setDependency(solver.diffusionEvent());
}


MutexSOSEvent::MutexSOSEvent(SOSSolver &solver,
                             string type,
                             string unit,
                             bool hasOutput,
                             bool storeValue) :
    SOSEvent(solver, type, unit, hasOutput, storeValue),
    m_mutexSolver(solver)
{

}
