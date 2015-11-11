#include "confinedconstantconcentration.h"

#include "../../sossolver.h"


ConfinedConstantConcentration::ConfinedConstantConcentration(SOSSolver &solver) :
    ConstantConcentration(solver)
{

}

void ConfinedConstantConcentration::registerHeightChange(const uint x,
                                                         const uint y,
                                                         const int value,
                                                         std::vector<DissolutionDeposition *> &affectedSurfaceReactions,
                                                         const uint nAffectedSurfaceReactions)
{
    (void) x;
    (void) y;
    (void) affectedSurfaceReactions;
    (void) nAffectedSurfaceReactions;

    m_deltaSum += value;
}

void ConfinedConstantConcentration::setupInitialConditions()
{
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
