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
    m_h0 = solver().averageHeight();
}

void GrowthSpeed::execute()
{
    const double &h = solver().averageHeight();

    if (cycle() != 0)
    {
        setValue((h - m_h0)/(solver().currentTime() - m_T0));
    }
}


void HeightRMS::execute()
{
    const double &mh = solver().averageHeight();

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


AutoCorrelationHeight::AutoCorrelationHeight(const SOSSolver &solver,
                                             const uint totalsamples,
                                             const uint interval) :
    SOSEvent(solver, "ach", "%",true),
    Observer(),
    m_xSpan(solver.length()/2),
    m_ySpan(solver.width()/2),
    m_autoCorrelationQuadrant(m_xSpan+1, m_ySpan+1),
    m_autoCorrelationSubQuadrant(m_xSpan == 0 ? 0 : m_xSpan,
                                 m_ySpan == 0 ? 0 : m_ySpan),
    m_totalsamples(totalsamples),
    m_interval(interval)
{
    if (m_ySpan > solver.width() || m_xSpan > solver.length())
    {
        throw std::logic_error( "invalid maxLength." );
    }
}

mat AutoCorrelationHeight::autoCorrelation() const
{
    mat autocorrelation(2*m_xSpan + 1, 2*m_ySpan + 1, fill::zeros);

    uint weight = solver().area()*m_count;

    //implement symmetry for 11 10 01 and -1-1 -10 0-1 steps
    for (uint dx = 0; dx <= m_xSpan; ++dx)
    {
        for (uint dy = 0; dy <= m_ySpan; ++dy)
        {
            autocorrelation(m_xSpan - (int)dx, m_ySpan - (int)dy)
                    = autocorrelation(m_xSpan + dx, m_ySpan + dy)
                    = m_autoCorrelationQuadrant(dx, dy)/weight;
        }
    }

    //implement symmetry for 1-1 and -11
    for (uint dx = 1; dx <= m_xSpan; ++dx)
    {
        for (uint dy = 1; dy <= m_ySpan; ++dy)
        {
            autocorrelation(m_xSpan + dx, m_ySpan - (int)dy)
                    = autocorrelation(m_xSpan - (int)dx, m_ySpan + dy)
                    = m_autoCorrelationSubQuadrant(dx - 1, dy - 1)/weight;
        }
    }

    return autocorrelation;

}

void AutoCorrelationHeight::updateFunction()
{
    int xn;
    int yn;

    const double &hMean = solver().averageHeight();

    //First we loop over all x,y that will create unique correlation pairs
    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            const int &h = solver().height(x, y);

            //correlate the point (x, y) with the rest of the system

            //positive directions
            for (uint dx = 0; dx <= m_xSpan; ++dx)
            {
                for (uint dy = 0; dy <= m_ySpan; ++dy)
                {
                    solver().boundaryLatticeTransform(xn, yn, x + dx, y + dy, h);

                    BADAssBool(!solver().isOutsideBox(xn, yn));

                    const int &hn = solver().height(xn, yn);

                    m_autoCorrelationQuadrant(dx, dy) += (hn-hMean)*(h-hMean);
                }
            }

            //diagonal downwards
            for (uint dx = 1; dx <= m_xSpan; ++dx)
            {
                for (uint dy = 1; dy <= m_ySpan; ++dy)
                {
                    solver().boundaryLatticeTransform(xn, yn, x + dx, y - (int)dy, h);

                    BADAssBool(!solver().isOutsideBox(xn, yn));

                    const int &hn = solver().height(xn, yn);

                    m_autoCorrelationSubQuadrant(dx-1, dy-1) += (hn-hMean)*(h-hMean);
                }
            }

        }
    }

    m_updateCount++;
}

void AutoCorrelationHeight::execute()
{
    setValue(m_updateCount/double(m_totalsamples)*100);

    if (cycle() % m_interval == 0)
    {
#ifdef DUMPCORR
        mat autocorr = autoCorrelation();
        autocorr.save("/tmp/autocorr.arma");
#endif
    }

    if (m_updateCount == m_totalsamples)
    {
        terminateLoop("AutoCorrelation Converged.");
    }
}

void AutoCorrelationHeight::initializeObserver(const Subjects &subject)
{
    (void) subject;

    m_autoCorrelationQuadrant.zeros();
    m_autoCorrelationSubQuadrant.zeros();
    m_count = 0;
    m_updateCount = 0;
}

void AutoCorrelationHeight::notifyObserver(const Subjects &subject)
{
    (void) subject;

    m_count++;

    if (((m_count-1) % m_interval) == 0)
    {
        updateFunction();
    }
}



void ConcentrationTracker::execute()
{
    setValue(solver().diffusionEvent().concentration());
}


AutoCorrelationHeightBF::AutoCorrelationHeightBF(const SOSSolver &solver,
                                                 const uint xSpan,
                                                 const uint ySpan,
                                                 const uint interval) :
    SOSEvent(solver, "autocorr", "", true),
    m_xSpan(xSpan == 0 ? solver.length()/2 : xSpan),
    m_ySpan(ySpan == 0 ? solver.width()/2 : ySpan),
    m_interval(interval),
    m_autoCorrelation(2*m_xSpan+1, 2*m_ySpan+1)
{

}

mat AutoCorrelationHeightBF::autoCorrelation() const
{
    return m_autoCorrelation/(m_count*solver().area());
}

void AutoCorrelationHeightBF::execute()
{
    setValue(m_count);

    if (cycle() % m_interval != 0)
    {
        return;
    }

    const double &hm = solver().averageHeight();

    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            const int &h = solver().height(x, y);

            for (int dx = -int(m_xSpan); dx <= int(m_xSpan); ++dx)
            {
                for (int dy = -int(m_ySpan); dy <= int(m_ySpan); ++dy)
                {
                    int xt, yt;
                    solver().boundaryLatticeTransform(xt, yt, int(x)+dx, int(y)+dy, h);

                    const int &ht = solver().height(xt, yt);

                    m_autoCorrelation(m_xSpan + dx, m_ySpan + dy) += (h-hm)*(ht-hm);
                }
            }
        }
    }

    m_count++;

    //#ifndef NDEBUG
    mat autocorr = autoCorrelation();
    autocorr.save("/tmp/autocorr" + to_string(m_count-1) + ".arma");
    //#endif
}

void AutoCorrelationHeightBF::initialize()
{
    m_autoCorrelation.zeros();
    m_count = 0;
}


AutoCorrelation1D::AutoCorrelation1D(const SOSSolver &solver) :
    SOSEvent(solver),
    m_acf(solver.length()/2)
{

}

vec AutoCorrelation1D::autoCorrelation() const
{
    return m_acf/((cycle()+1)*solver().length()/2);
}

void AutoCorrelation1D::execute()
{
    const double &mh = solver().averageHeight();

    for (uint dx = 0; dx < solver().length()/2; ++dx)
    {
        for (uint x = 0; x < solver().length()/2; ++x)
        {
            m_acf(dx) += (solver().height(x, 0) - mh)*(solver().height(x+dx,0) - mh);
        }
    }

    m_acf.save("/tmp/acf.arma");
}

void AutoCorrelation1D::initialize()
{
    m_acf.zeros();
}

void FreeVolumeTracker::execute()
{
    setValue(solver().freeVolume());
}
