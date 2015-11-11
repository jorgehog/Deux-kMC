#include "constantconfinement.h"

#include "../../sossolver.h"


ConstantConfinement::ConstantConfinement(SOSSolver &solver, const double height) :
    ConfiningSurface(solver, "ConstantConfinement"),
    FixedSurface(solver, height)
{

}

ConstantConfinement::~ConstantConfinement()
{

}



void ConstantConfinement::registerHeightChange(const uint x, const uint y, const int value, std::vector<DissolutionDeposition *> &affectedSurfaceReactions, const uint nAffectedSurfaceReactions)
{
    (void) x;
    (void) y;
    (void) value;
    (void) affectedSurfaceReactions;
    (void) nAffectedSurfaceReactions;

    double newHeight = solver().averageHeight() + m_h;

    if (solver().heights().max() <= newHeight - 1)
    {
        setHeight(newHeight);
    }
}


void ConstantConfinement::initialize()
{
    m_h = height() - solver().averageHeight();
}

