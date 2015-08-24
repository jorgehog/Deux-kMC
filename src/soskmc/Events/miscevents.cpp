#include "miscevents.h"
#include "sosreaction.h"
#include "../kmcsolver/boundary/boundary.h"

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

double SurfaceSize::relativeHeightSum(const uint x, const uint y) const
{
    const uint right = solver().rightSite(x);
    const uint top = solver().topSite(y);

    if (solver().isOutsideBoxSingle(right, 0) || solver().isOutsideBoxSingle(top, 1))
    {
        return 0;
    }

    double xDirection = abs(solver().height(x, y) - solver().height(right, y));
    double yDirection = abs(solver().height(x, y) - solver().height(x, top));

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
            value += relativeHeightSum(x, y);
        }
    }

    value /= (solver().area());

    return value;

}

void SurfaceSize::execute()
{
    m_sum += solver().currentTimeStep()*getLocalValue();

    setValue(getLocalValue());
}

void SurfaceSize::reset()
{
    for (const auto &xy : solver().changedSurfaceSites())
    {
        const uint &x = xy.first;
        const uint &y = xy.second;

        const int left = solver().leftSite(x);
        const int bottom = solver().bottomSite(y);

        updateRelativeHeight(x, y);

        if (!solver().boundary(0, 0)->isBlocked(left))
        {
            updateRelativeHeight(left, y);
        }

        if (!solver().boundary(1, 0)->isBlocked(bottom))
        {
            updateRelativeHeight(x, bottom);
        }

    }

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


void AverageHeight::initialize()
{
    setValue(getValue());
}

void AverageHeight::execute()
{
    setValue(getValue());
}

double AverageHeight::getValue() const
{
    return accu(solver().heights())/(double)solver().area();
}


void NNeighbors::updateNNeighbors(const uint x, const uint y)
{
    uint newValue = solver().nNeighbors(x, y);

    m_localValue += ((int)newValue - (int)m_nNeighbors(x, y));

    m_nNeighbors(x, y) = newValue;
}

double NNeighbors::bruteForceValue() const
{
    double value = 0;

    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            value += solver().nNeighbors(x, y);
        }
    }

    return value/solver().area();
}

void NNeighbors::initialize()
{
    m_localValue = 0;

    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            uint n = solver().nNeighbors(x, y);

            m_nNeighbors(x, y) = n;
            m_localValue += n;
        }
    }
}

void NNeighbors::execute()
{
    setValue(m_localValue/solver().area());
}

void NNeighbors::reset()
{

    for (const auto &xy : solver().changedSurfaceSites())
    {
        const uint &x = xy.first;
        const uint &y = xy.second;

        //no change
        if (m_nNeighbors(x, y) == solver().nNeighbors(x, y))
        {
            continue;
        }

        const uint left = solver().leftSite(x);
        const uint right = solver().rightSite(x);
        const uint bottom = solver().bottomSite(y);
        const uint top = solver().topSite(y);

        updateNNeighbors(x, y);

        if (!solver().isOutsideBoxSingle(left, 0))
        {
            updateNNeighbors(left, y);
        }

        if (!solver().isOutsideBoxSingle(right, 0))
        {
            updateNNeighbors(right, y);
        }

        if (!solver().isOutsideBoxSingle(bottom, 1))
        {
            updateNNeighbors(x, bottom);
        }

        if (!solver().isOutsideBoxSingle(top, 1))
        {
            updateNNeighbors(x, top);
        }

    }
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

        if (reaction->isAllowed())
        {
            BADAss(reaction->rate(), !=, 0);

            BADAssClose(reaction->rate(), reaction->rateExpression(), 1E-5, "error in rate updating.", [&] ()
            {
                BADAssSimpleDump(reaction->rate(), reaction->rateExpression());
            });
        }

    }
}


void DumpHeightSlice::initialize()
{
    if (m_axis == 0)
    {
        m_heights.set_size(solver().heights().n_rows);
    }

    else
    {
        m_heights.set_size(solver().heights().n_cols);
    }
}

void DumpHeightSlice::execute()
{
    if (cycle() % m_nCyclesPerOutput != 0)
    {
        return;
    }

    if (m_axis == 0)
    {
        m_heights = solver().heights().col(m_slicePosition);
    }

    else
    {
        m_heights = solver().heights().row(m_slicePosition);
    }

    m_heights.save(m_path + "/heights.arma");
}


void SurfaceVariance::execute()
{
    m_s2 += solver().currentTimeStep()*pow(dependency<SurfaceSize>("SurfaceSize")->getLocalValue(), 2);

    if (cycle() != 0)
    {
        setValue(sqrt(m_s2/(solver().currentTime()-m_T0) - pow(dependency<SurfaceSize>("SurfaceSize")->timeAverage(), 2)));
    }
}



void GrowthSpeed::initialize()
{
    m_T0 = solver().currentTime();
    m_h0 = dependency<AverageHeight>("AverageHeight")->value();
}

void GrowthSpeed::execute()
{
    const double &h = dependency<AverageHeight>("AverageHeight")->value();

    if (cycle() != 0)
    {
        setValue((h - m_h0)/(solver().currentTime() - m_T0));
    }
}


void HeightRMS::execute()
{
    const double &mh = dependency<AverageHeight>("AverageHeight")->value();

    double rms = 0;

    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            double deviation = solver().height(x, y) - mh;

            rms += deviation*deviation;
        }
    }

    rms /= solver().area();

    setValue(sqrt(rms));
}
