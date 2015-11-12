#include "constantconcentration.h"

ConstantConcentration::ConstantConcentration(SOSSolver &solver) :
    ConcentrationProfile(solver, [] (const uint x, const uint y) {(void) x; (void) y; return 1.0;})
{

}

ConstantConcentration::~ConstantConcentration()
{

}
