#include "confinedconstantconcentration.h"

#include "../../sossolver.h"
#include "../confiningsurface/confiningsurface.h"

ConfinedConstantConcentration::ConfinedConstantConcentration(SOSSolver &solver) :
    ConstantConcentration(solver),
    m_fix(false)
{

}

void ConfinedConstantConcentration::notifyObserver(const Subjects &subject)
{
    if (subject == Subjects::SOLVER)
    {
        if (solver().currentSurfaceChange().type == ChangeTypes::Single)
        {
            m_deltaSum += solver().currentSurfaceChange().value;
        }
    }
}

void ConfinedConstantConcentration::initializeObserver(const Subjects &subject)
{
    (void) subject;

    m_N0 = solver().concentration()*solver().freeVolume();

    m_deltaSum = 0;
}

void ConfinedConstantConcentration::reset()
{
    if (m_fix)
    {
        return;
    }

    const double Vf = solver().freeVolume();

    if (m_deltaSum < m_N0 && Vf != 0)
    {
        solver().setConcentration((m_N0 - m_deltaSum)/Vf);
    }

    else
    {
        solver().setZeroConcentration();
    }
}
