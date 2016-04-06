#include "confinedconstantconcentration.h"

#include "../../sossolver.h"
#include "../confiningsurface/confiningsurface.h"

ConfinedConstantConcentration::ConfinedConstantConcentration(SOSSolver &solver) :
    ConstantConcentration(solver)
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

    else
    {
        const ConfiningSurface &cse = solver().confiningSurfaceEvent();

        const double dh = cse.height() - cse.currentConfinementChange().prevHeight;

        m_V0 += solver().area()*dh;
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
    BADAssClose(m_V0 - m_deltaSum, solver().freeVolume(), 1E-10);

    if (m_deltaSum < m_N0 && m_deltaSum < m_V0)
    {
        solver().setConcentration((m_N0 - m_deltaSum)/(m_V0 - m_deltaSum));
    }

    else
    {
        solver().setZeroConcentration();
    }
}
