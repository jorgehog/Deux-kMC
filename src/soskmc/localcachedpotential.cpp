#include "localcachedpotential.h"

#include "sossolver.h"


LocalCachedPotential::LocalCachedPotential(SOSSolver &solver) :
    LocalPotential(solver),
    m_potentialValues(solver.length(), solver.width())
{

}

double LocalCachedPotential::sum() const
{
    return arma::accu(m_potentialValues);
}

double LocalCachedPotential::potential(const uint x, const uint y) const
{
    BADAssClose(m_potentialValues(x, y), potentialFunction(x, y), 1E-5, "incorrect cache", [&] ()
    {
        BADAssSimpleDump(m_potentialValues(x, y),
                         potentialFunction(x, y));
    });

    return m_potentialValues(x, y);
}

void LocalCachedPotential::initializeObserver(const Subjects &subject)
{
    (void) subject;

    //this will be called for every subject... not ideal, but fine since
    //it is called only once per subject, unlike notify
    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            m_potentialValues(x, y) = potentialFunction(x, y);
        }
    }
}
