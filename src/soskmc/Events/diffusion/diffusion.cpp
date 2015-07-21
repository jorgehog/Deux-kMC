#include "diffusion.h"

Diffusion::Diffusion(SOSSolver &solver,
                     string type,
                     string unit,
                     bool hasOutput,
                     bool storeValue) :
    SolidOnSolidEvent(solver, type, unit, hasOutput, storeValue)
{
    solver.setDiffusionEvent(*this);
}

Diffusion::~Diffusion()
{

}


