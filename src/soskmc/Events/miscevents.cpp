#include "miscevents.h"
#include "sosreaction.h"
#include "../kmcsolver/boundary/boundary.h"
#include "diffusion/diffusion.h"

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
    const int &h = solver().height(x, y);

    const uint right = solver().rightSite(x, y, h);
    const uint top = solver().topSite(x, y, h);

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

        const int &h = solver().height(x, y);

        const int left = solver().leftSite(x, y, h);
        const int bottom = solver().bottomSite(x, y, h);

        updateRelativeHeight(x, y);

        if (!solver().isOutsideBoxSingle(left, 0))
        {
            updateRelativeHeight(left, y);
        }

        if (!solver().isOutsideBoxSingle(bottom, 1))
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
    m_h0 = solver().averageHeight();
}

void AverageHeight::execute()
{
    setValue(getValue());
}

double AverageHeight::getValue() const
{
    BADAssClose(solver().averageHeight(), accu(solver().heights())/double(solver().area()), 1E-3);
    return solver().averageHeight() - m_h0;
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
        const int &h = solver().height(x, y);

        //no change
        if (m_nNeighbors(x, y) == solver().nNeighbors(x, y))
        {
            continue;
        }

        const uint left = solver().leftSite(x, y, h);
        const uint right = solver().rightSite(x, y, h);
        const uint bottom = solver().bottomSite(x, y, h);
        const uint top = solver().topSite(x, y, h);

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

void RateChecker::execute()
{
    for (uint n = 0; n < m_solver.numberOfReactions(); ++n)
    {
        Reaction *reaction = m_solver.getReaction(n);

        if (reaction->isAllowed())
        {
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


AutoCorrelationHeight::AutoCorrelationHeight(const SOSSolver &solver, const uint xSpan, const uint ySpan) :
    SOSEvent(solver),
    m_xSpan(xSpan == 0 ? solver.length()/2 : xSpan),
    m_ySpan(ySpan == 0 ? solver.width()/2 : ySpan),
    m_autoCorrelationQuadrant(m_xSpan, m_ySpan),
    m_autoCorrelationSubQuadrant(m_xSpan == 0 ? 0 : m_xSpan - 1,
                                 m_ySpan == 0 ? 0 : m_ySpan - 1)
{
    if (m_ySpan > solver.width() || m_xSpan > solver.length())
    {
        throw std::logic_error( "invalid maxLength." );
    }
}

mat AutoCorrelationHeight::autoCorrelation() const
{
    mat autocorrelation(2*m_xSpan - 1, 2*m_ySpan - 1, fill::zeros);

    uint weight = m_xSpan*m_ySpan*(cycle() + 1);

    //implement symmetry for 11 10 01 and -1-1 -10 0-1 steps
    for (uint dx = 0; dx < m_xSpan; ++dx)
    {
        for (uint dy = 0; dy < m_ySpan; ++dy)
        {
            autocorrelation(m_xSpan - dx - 1, m_ySpan - dy - 1)
                    = autocorrelation(m_xSpan + dx - 1, m_ySpan + dy - 1)
                    = m_autoCorrelationQuadrant(dx, dy)/weight;
        }
    }

    uint subweight = (m_xSpan - 1)*(m_ySpan - 1)*(cycle() + 1);

    //implement symmetry for 1-1 and -11
    for (uint dx = 1; dx < m_xSpan; ++dx)
    {
        for (uint dy = 1; dy < m_ySpan; ++dy)
        {
            autocorrelation(m_xSpan + dx - 1, m_ySpan - dy - 1)
                    = autocorrelation(m_xSpan - dx - 1, m_ySpan + dy - 1)
                    = m_autoCorrelationSubQuadrant(dx - 1, dy - 1)/subweight;
        }
    }

    return autocorrelation;

}

void AutoCorrelationHeight::execute()
{
    int xn;
    int yn;

    const double &hMean = solver().averageHeight();

    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            const int &h = solver().height(x, y);

            //correlate the point (x, y) with the rest of the system

            //positive directions
            for (uint dx = 0; dx < m_xSpan; ++dx)
            {
                for (uint dy = 0; dy < m_ySpan; ++dy)
                {
                    solver().boundaryLatticeTransform(xn, yn, x + dx, y + dy, h);

                    if (!solver().isOutsideBox(xn, yn))
                    {
                        const int &hn = solver().height(xn, yn);

                        m_autoCorrelationQuadrant(dx, dy) += (hn-hMean)*(h-hMean);

                        //derp: if something is outside box, weighting is no longer area*cycles
                    }

                }
            }

            //diagonal downwards
            for (uint dx = 1; dx < m_xSpan; ++dx)
            {
                for (uint dy = 1; dy < m_ySpan; ++dy)
                {
                    solver().boundaryLatticeTransform(xn, yn, x + dx, y - (int)dy, h);

                    if (!solver().isOutsideBox(xn, yn))
                    {
                        const int &hn = solver().height(xn, yn);

                        m_autoCorrelationSubQuadrant(dx-1, dy-1) += (hn-hMean)*(h-hMean);

                        //derp: if something is outside box, weighting is no longer area*cycles
                    }
                }
            }

        }
    }

    //#ifndef NDEBUG
    mat autocorr = autoCorrelation();
    autocorr.save("/tmp/autocorr.arma");
    //#endif
}

void AutoCorrelationHeight::initialize()
{
    m_autoCorrelationQuadrant.zeros();
    m_autoCorrelationSubQuadrant.zeros();
}


void ConcentrationTracker::execute()
{
    setValue(solver().diffusionEvent().concentration());
}
