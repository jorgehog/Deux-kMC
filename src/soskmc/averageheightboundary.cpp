#include "averageheightboundary.h"
#include "dissolutiondeposition.h"


AverageHeightBoundary::AverageHeightBoundary(SOSSolver &solver,
                                             const uint cutoff,
                                             const uint dim,
                                             const uint span,
                                             const uint yspan,
                                             Boundary::orientations orientation,
                                             const double location) :
    SOSBoundary(solver, orientation),
    m_mutexSolver(solver),
    m_cutoff(cutoff == 0 ? span : cutoff),
    m_dim(dim),
    m_span(span),
    m_yspan(yspan),
    m_location(location)
{
    getStartsAndEnds(m_x0, m_y0, m_x1, m_y1);

    m_mutexSolver.registerHeightConnecter(this);
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
            m_mutexSolver.registerChangedSite(m_location, y);
        }

        else
        {
            m_mutexSolver.registerChangedSite(y, m_location);
        }
    }
}

void AverageHeightBoundary::affectSolutionSites(const int z)
{
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
        if (m_orientation == orientations::FIRST)
        {
            x0 = 0;
            x1 = m_cutoff;
        }

        else
        {
            x0 = m_span - m_cutoff;
            x1 = m_span;
        }

        y0 = 0;
        y1 = m_yspan;
    }

    else
    {
        if (m_orientation == orientations::FIRST)
        {
            y0 = 0;
            y1 = m_cutoff;
        }

        else
        {
            y0 = m_span - m_cutoff;
            y1 = m_span;
        }

        x0 = 0;
        x1 = m_yspan;
    }
}

double AverageHeightBoundary::transformCoordinate(const double xi, const double xj, const double xk) const
{
    (void) xj;
    (void) xk;

    return xi;
}

bool AverageHeightBoundary::isBlocked(const double xi, const double xj, const double xk) const
{
    (void) xj;

    const double &averageHeight = average();

    if (m_orientation == orientations::FIRST)
    {
        return (xi < m_location) && (xk <= averageHeight);
    }

    else
    {
        return (xi > m_location) && (xk <= averageHeight);
    }
}

std::vector<double> AverageHeightBoundary::imagesOf(const double xi, const double xj, const double xk) const
{
    (void) xi;
    (void) xj;
    (void) xk;

    return {};
}


void AverageHeightBoundary::registerHeightChange(const uint x, const uint y, const int value, std::vector<DissolutionDeposition *> &affectedSurfaceReactions, const uint nAffectedSurfaceReactions)
{
    (void) affectedSurfaceReactions;
    (void) nAffectedSurfaceReactions;

    if (isInsideCutoff(x, y))
    {
        m_average += double(value)/cutoffArea();
    }

    BADAssClose(m_average, calcAverage(), 1E-3);

    int height = round(m_average);
    int prevHeight = round(m_prevAverage);

    m_prevAverage = m_average;

    updateSites(height, prevHeight);
}

void AverageHeightBoundary::setupInitialConditions()
{
    m_average = calcAverage();

    double prevPrev = m_prevAverage;
    m_prevAverage = m_average;

    updateSites(round(m_prevAverage), round(prevPrev));
}
