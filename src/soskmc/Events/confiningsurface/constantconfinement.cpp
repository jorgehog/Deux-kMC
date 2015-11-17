#include "constantconfinement.h"

#include "../../sossolver.h"


ConstantConfinement::ConstantConfinement(SOSSolver &solver, const double height) :
    ConfiningSurface(solver, "ConstantConfinement", "", true, true),
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

    if (maxHeight <= newHeight - 1)
    {
//        cout << "changing height "<< newHeight << endl;
        setHeight(newHeight);
//        cout << "height changed." << endl;
//        cout << "it crashes because the surface has been shifted down before the particle has been released, so it is not 'cleaned up' by the notification from the confined surface." << endl;
    }

//    else
//    {
//        cout << "need to dissolve!" << setprecision(16) << fixed << height() << endl;
//        BADAss(maxHeight, <=, height() - 1);
//        mutexSolver().registerHeightChange(xi, yi, -1);
//    }
}


void ConstantConfinement::initialize()
{
    m_h = height() - solver().averageHeight();
}

void ConstantConfinement::execute()
{
    setValue((height() - solver().averageHeight() - m_h));
}

