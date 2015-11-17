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



void ConstantConfinement::notifyObserver()
{
    double newHeight = solver().averageHeight() + m_h;

    uint xi, yi;
    int maxHeight = solver().heights().max(xi, yi);

    if (maxHeight <= newHeight - 1)
    {
        setHeight(newHeight);
    }

    else
    {
        mutexSolver().registerHeightChange(xi, yi, -1);
    }
}


void ConstantConfinement::initialize()
{
    m_h = height() - solver().averageHeight();
}

