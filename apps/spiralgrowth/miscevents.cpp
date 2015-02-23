#include "miscevents.h"



void SurfaceSize::execute()
{
    m_localValue = 0;

    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            m_localValue += 0.5*abs(solver().height(x, y) - solver().height(solver().rightSite(x), y));
            m_localValue += 0.5*abs(solver().height(x, y) - solver().height(x, solver().topSite(y)));
        }
    }

    m_localValue /= (solver().length()*solver().width());

    m_sum += solver().currentTimeStep()*m_localValue;

    setValue(m_sum/(solver().currentTime() - m_T0));
}


void DumpHeights3D::initialize()
{
    m_L = solver().length();
    m_W = solver().width();
    m_N = m_L*m_W;
}

void DumpHeights3D::execute()
{
    if (cycle() % 100 != 0)
    {
        return;
    }

    const imat &heights = solver().heights();

    const int max = heights.max();
    const int min = heights.min();
    const uint maxSpan = max - min;

    m_writer.setSystemSize(m_L, m_W, maxSpan, 0, 0, 0);

    m_writer.initializeNewFile(cycle()/100);

    for (uint x = 0; x < m_L; ++x)
    {
        for (uint y = 0; y < m_W; ++y)
        {
            uint zSpan = heights(x, y) - min;

            for (uint z = 0; z < zSpan; ++z)
            {
                m_writer << x
                         << y
                         << z
                         << 0;
            }

            m_writer << x
                     << y
                     << zSpan
                     << solver().nNeighbors(x, y);
        }
    }

    m_writer.finalize();

}


void AverageHeight::execute()
{
    setValue(accu(solver().heights())/(double)solver().area());
}


void NNeighbors::execute()
{
}




RateChecker::RateChecker(const KMCSolver &solver) :
    LatticeEvent("rateChecker"),
    m_solver(solver)
{

}

void RateChecker::reset()
{
    for (uint n = 0; n < m_solver.numberOfReactions(); ++n)
    {
        Reaction *reaction = m_solver.getReaction(n);

        if (!reaction->isAllowed())
        {
            BADAss(reaction->rate(), ==, 0);
        }

        else
        {
            BADAss(reaction->rate(), !=, 0);

            BADAssClose(reaction->rate(), reaction->rateExpression(), 1E-5, "error in rate updating.", [&] ()
            {
                BADAssSimpleDump(reaction->rate(), reaction->rateExpression());
            });
        }

    }
}
