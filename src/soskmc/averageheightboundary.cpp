#include "averageheightboundary.h"
#include "dissolutiondeposition.h"
#include "Events/diffusion/diffusion.h"


AverageHeightBoundary::AverageHeightBoundary(SOSSolver &solver,
                                             const uint averageHeightDepth,
                                             const uint dim,
                                             const uint span,
                                             const uint yspan,
                                             Boundary::orientations orientation,
                                             const uint location) :
    SOSBoundary(solver, orientation),
    m_averageHeightDepth(averageHeightDepth == 0 ? span : averageHeightDepth),
    m_dim(dim),
    m_span(span),
    m_yspan(yspan),
    m_location(location)
{
    getStartsAndEnds(m_x0, m_y0, m_x1, m_y1);

    solver.registerObserver(this);
}

AverageHeightBoundary::~AverageHeightBoundary()
{

}

void AverageHeightBoundary::affectSurfaceSites()
{

    for (uint y = 0; y < m_yspan; ++y)
    {
        if (m_dim == 0)
        {
            solver().registerChangedSite(m_location, y);
        }

        else
        {
            solver().registerChangedSite(y, m_location);
        }
    }
}


double AverageHeightBoundary::calcAverage() const
{
    double average = 0;

    for (uint x = m_x0; x < m_x1; ++x)
    {
        for (uint y = m_y0; y < m_y1; ++y)
        {
            average += solver().height(x, y);
        }
    }

    return average / cutoffArea();
}

bool AverageHeightBoundary::isInsideCutoff(const int x, const int y) const
{
    if (m_dim == 0)
    {
        return (x >= (int)m_x0) && (x < (int)m_x1);
    }

    else
    {
        return (y >= (int)m_y0) && (y < (int)m_y1);
    }
}

void AverageHeightBoundary::getStartsAndEnds(uint &x0, uint &y0, uint &x1, uint &y1)
{
    if (m_dim == 0)
    {
        if (orientation() == orientations::FIRST)
        {
            x0 = 1;
            x1 = m_averageHeightDepth + 1;
        }

        else
        {
            x0 = m_span - m_averageHeightDepth - 1;
            x1 = m_span - 1;
        }

        y0 = 0;
        y1 = m_yspan;
    }

    else
    {
        if (orientation() == orientations::FIRST)
        {
            y0 = 1;
            y1 = m_averageHeightDepth + 1;
        }

        else
        {
            y0 = m_span - m_averageHeightDepth - 1;
            y1 = m_span - 1;
        }

        x0 = 0;
        x1 = m_yspan;
    }
}

double AverageHeightBoundary::transformContinousCoordinate(const double xi, const double xj, const double xk) const
{
    return transformFunction(xi, xj, xk);
}

int AverageHeightBoundary::transformLatticeCoordinate(const int xi, const int xj, const int xk) const
{
    return transformFunction(xi, xj, xk);
}

bool AverageHeightBoundary::isBlockedContinous(const double xi, const double xj, const double xk) const
{
    return blockedFunction(xi, xj, xk);
}

bool AverageHeightBoundary::isBlockedLattice(const int xi, const int xj, const int xk) const
{
    return blockedFunction(xi, xj, xk);
}


void AverageHeightBoundary::closestImage(const double xi, const double xj, const double xk,
                                         const double xti, const double xtj, const double xtk,
                                         double &dxi, double &dxj, double &dxk) const
{
    noImage(xi, xj, xk, xti, xtj, xtk, dxi, dxj, dxk);
}


void AverageHeightBoundary::notifyObserver(const Subjects &subject)
{
    (void) subject;

    const uint &x = solver().currentSurfaceChange().x;
    const uint &y = solver().currentSurfaceChange().y;
    const int &value = solver().currentSurfaceChange().value;

    if (solver().currentSurfaceChange().type == ChangeTypes::Single)
    {
        if (isInsideCutoff(x, y))
        {
            m_average += double(value)/cutoffArea();
        }
    }

    else
    {
        const int &x1 = solver().currentSurfaceChange().x1;
        const int &y1 = solver().currentSurfaceChange().y1;

        const bool startInside = isInsideCutoff(x, y);
        const bool endInside = isInsideCutoff(x1, y1);

        if (startInside && !endInside)
        {
            m_average -= 1./cutoffArea();
        }

        else if (!startInside && endInside)
        {
            m_average += 1./cutoffArea();
        }
    }

    BADAssClose(m_average, calcAverage(), 1E-3);

    int height = round(m_average);
    int prevHeight = round(m_prevAverage);

    m_prevAverage = m_average;

    if (height != prevHeight)
    {
        affectSurfaceSites();
    }
}

void AverageHeightBoundary::initializeObserver(const Subjects &subject)
{
    (void) subject;

    m_average = calcAverage();

    double prev;
    if (solver().cycle() == 0)
    {
        prev = m_average-1;
    }
    else
    {
        prev = m_prevAverage;
    }
    m_prevAverage = m_average;

    if (round(m_prevAverage) != round(prev))
    {
        affectSurfaceSites();
    }
}
