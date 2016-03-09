#include "averageheightboundary.h"
#include "dissolutiondeposition.h"
#include "Events/diffusion/diffusion.h"


AverageHeightBoundary::AverageHeightBoundary(SOSSolver &solver,
                                             const uint averageHeightDepth,
                                             const uint dim,
                                             const uint span,
                                             const uint yspan,
                                             Boundary::orientations orientation,
                                             const double location) :
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

void AverageHeightBoundary::affectSolutionSites(const int z)
{
    BADAssBool(!solver().diffusionEvent().hasDiscreteParticles());
    (void) z;
    //derp
}

void AverageHeightBoundary::updateSites(const int height, const int prevHeight)
{
    int delta;
    if (height == prevHeight)
    {
        return;
    }

    else if (height < prevHeight)
    {
        delta = -1;
    }

    else
    {
        delta = 1;
    }

    affectSurfaceSites();

    for (int z = prevHeight; z != height + delta; z += delta)
    {
        affectSolutionSites(z);
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

bool AverageHeightBoundary::isInsideCutoff(const uint x, const uint y) const
{
    if (m_dim == 0)
    {
        return (x >= m_x0) && (x < m_x1);
    }

    else
    {
        return (y >= m_y0) && (y < m_y1);
    }
}

void AverageHeightBoundary::getStartsAndEnds(uint &x0, uint &y0, uint &x1, uint &y1)
{
    if (m_dim == 0)
    {
        if (orientation() == orientations::FIRST)
        {
            x0 = 0;
            x1 = m_averageHeightDepth;
        }

        else
        {
            x0 = m_span - m_averageHeightDepth;
            x1 = m_span;
        }

        y0 = 0;
        y1 = m_yspan;
    }

    else
    {
        if (orientation() == orientations::FIRST)
        {
            y0 = 0;
            y1 = m_averageHeightDepth;
        }

        else
        {
            y0 = m_span - m_averageHeightDepth;
            y1 = m_span;
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
        const uint &x1 = solver().currentSurfaceChange().x1;
        const uint &y1 = solver().currentSurfaceChange().y1;

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

    updateSites(height, prevHeight);
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

    updateSites(round(m_prevAverage), round(prev));
}
