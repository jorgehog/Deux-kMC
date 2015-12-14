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



void ConstantConfinement::notifyObserver(const Subjects &subject)
{
    (void) subject;

    double newHeight = solver().averageHeight() + m_h;

    uint xi, yi;
    int maxHeight = solver().heights().max(xi, yi);

    BADAssEqual(maxHeight, solver().height(xi, yi));

    //there is no good solution if this criteria is violated.
    //we might rearrange the surface, but then we will decrease the mean height and
    //this function will be call repeatadly, and in some cases indefenetely.
    if (maxHeight <= newHeight - 1)
    {
        setHeight(newHeight);
    }
}


void ConstantConfinement::initialize()
{
    m_h = height() - solver().averageHeight();
}

void ConstantConfinement::execute()
{
    setValue((height() - solver().averageHeight() - m_h));
}

