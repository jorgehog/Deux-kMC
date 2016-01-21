#include "localcachedpotential.h"

#include "sossolver.h"


LocalCachedPotential::LocalCachedPotential(SOSSolver &solver) :
    LocalPotential(solver),
    m_potentialValues(solver.length(), solver.width())
{

}

double LocalCachedPotential::potential(const uint x, const uint y) const
{
    BADAssClose(m_potentialValues(x, y), potentialFunction(x, y), 1E-5);
    return m_potentialValues(x, y);
}

void LocalCachedPotential::initializeObserver(const Subjects &subject)
{
    if (subject == Subjects::SOLVER)
    {
        for (uint x = 0; x < solver().length(); ++x)
        {
            for (uint y = 0; y < solver().width(); ++y)
            {
                m_potentialValues(x, y) = potentialFunction(x, y);
            }
        }
    }
}
