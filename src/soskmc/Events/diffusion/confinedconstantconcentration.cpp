#include "confinedconstantconcentration.h"

#include "../../sossolver.h"


ConfinedConstantConcentration::ConfinedConstantConcentration(SOSSolver &solver) :
    ConstantConcentration(solver)
{

}

void ConfinedConstantConcentration::notifyObserver(const Subjects &subject)
{
    (void) subject;

    if (solver().currentSurfaceChange().type == ChangeTypes::Single)
    {
        m_deltaSum += solver().currentSurfaceChange().value;
    }
}

void ConfinedConstantConcentration::initializeObserver(const Subjects &subject)
{
    (void) subject;

    m_V0 = solver().freeVolume();
    m_N0 = solver().concentration()*m_V0;

    m_deltaSum = 0;
}

void ConfinedConstantConcentration::reset()
{
    if (m_deltaSum < m_N0 && m_deltaSum < m_V0)
    {
        solver().setConcentration((m_N0 - m_deltaSum)/(m_V0 - m_deltaSum));
    }

    else
    {
        solver().setZeroConcentration();
    }
}
