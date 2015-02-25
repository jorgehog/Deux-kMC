#include "miscevents.h"
#include "solidonsolidreaction.h"


void SurfaceSize::initialize()
{
    m_sum = 0;
    m_T0 = solver().currentTime() - solver().currentTimeStep();

    m_localValue = 0;

    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            double _relativeHeightSum = relativeHeightSum(x, y);

            m_localValue += _relativeHeightSum;
            m_relativeHeightSums(x, y) = _relativeHeightSum;
        }
    }
}

double SurfaceSize::relativeHeightSum(const uint x, const uint y)
{
    double xDirection = abs(solver().height(x, y) - solver().height(solver().rightSite(x), y));
    double yDirection = abs(solver().height(x, y) - solver().height(x, solver().topSite(y)));

    return 0.5*(xDirection + yDirection);
}

void SurfaceSize::updateRelativeHeight(const uint x, const uint y)
{
    double newValue = relativeHeightSum(x, y);

    m_localValue += (newValue - m_relativeHeightSums(x, y));

    m_relativeHeightSums(x, y) = newValue;
}

double SurfaceSize::bruteForceValue() const
{
    double value = 0;

    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            value += 0.5*abs(solver().height(x, y) - solver().height(solver().rightSite(x), y));
            value += 0.5*abs(solver().height(x, y) - solver().height(x, solver().topSite(y)));
        }
    }

    value /= (solver().area());

    return value;

}

void SurfaceSize::execute()
{
    m_sum += solver().currentTimeStep()*getLocalValue();

    setValue(m_sum/(solver().currentTime() - m_T0));
}

void SurfaceSize::reset()
{
    const DiffusionDeposition *currentReaction = dynamic_cast<const DiffusionDeposition*>(solver().selectedReaction());

    const uint &x = currentReaction->x();
    const uint &y = currentReaction->y();

    const uint leftSite = solver().leftSite(x);
    const uint bottomSite = solver().bottomSite(y);

    updateRelativeHeight(x, y);
    updateRelativeHeight(leftSite, y);
    updateRelativeHeight(x, bottomSite);

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


void NNeighbors::initialize()
{
    m_sum = 0;
}

void NNeighbors::execute()
{
    m_localValue = 0;

    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            m_localValue += solver().nNeighbors(x, y);
        }
    }

    m_localValue /= solver().area();

    m_sum += m_localValue;

    setValue(value());
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
