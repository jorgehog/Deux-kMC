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

    if (solver().currentSurfaceChange().type == ChangeTypes::Double)
    {
        return;
    }

    const double &avgh = solver().averageHeight();
    const double newHeight = avgh + m_h;

    uint xi, yi;
    const int maxHeight = solver().heights().max(xi, yi);

    BADAssEqual(maxHeight, solver().height(xi, yi));

    //there is no good solution if this criteria is violated.
    //we might rearrange the surface, but then we will decrease the mean height and
    //this function will be call repeatadly, and in some cases indefenetely.
    if (maxHeight <= newHeight - 1)
    {
        setHeight(newHeight);
    }

    else
    {
        setHeight(maxHeight + 1);
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

