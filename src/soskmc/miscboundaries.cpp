#include "miscboundaries.h"

#include "sossolver.h"
#include "dissolutiondeposition.h"

BoundaryTrackingDevice::BoundaryTrackingDevice(SOSSolver &solver,
                                               const uint location,
                                               const uint dim,
                                               const uint depth,
                                               const bool affectOnChange) :
    Observer(),
    m_solver(solver),
    m_location(location),
    m_dim(dim),
    m_depth(depth),
    m_area(solver.orthExtent(dim)*depth),
    m_affectOnChange(affectOnChange)
{
    if (location == 0)
    {
        m_x0 = 1;
    }

    else
    {
        m_x0 = solver.extent(dim) - depth - 1;
    }

    m_x1 = m_x0 + depth;
}

bool BoundaryTrackingDevice::pointIsInsideArea(const int x, const int y) const
{
    if (solver().isOutsideBox(x, y))
    {
        return false;
    }

    uint ux;
    if (dim() == 0)
    {
        ux = uint(x);
    }

    else
    {
        ux = uint(y);
    }

    return (ux >= x0() && ux < x1());
}

void BoundaryTrackingDevice::onAverageChange()
{
    if (!m_affectOnChange)
    {
        return;
    }

    SurfaceReaction *r;
    for (uint boundarySite = 0; boundarySite < solver().orthExtent(dim()); ++boundarySite)
    {
        if (dim() == 0)
        {
            r = &solver().surfaceReaction(m_location, boundarySite);
        }

        else
        {
            r = &solver().surfaceReaction(boundarySite, m_location);
        }

        solver().registerAffectedReaction(r);
    }
}


TrackLineAverage::TrackLineAverage(SOSSolver &solver,
                                   const uint location,
                                   const uint dim,
                                   const uint depth,
                                   const bool affectOnChange) :
    BoundaryTrackingDevice(solver, location, dim, depth, affectOnChange),
    m_averages(solver.orthExtent(dim))
{

}

double TrackLineAverage::bruteForceAverage(const uint boundarySite) const
{
    int hAvg = 0;
    for (uint lineSite = x0(); lineSite < x1(); ++lineSite)
    {
        if (dim() == 0)
        {
            hAvg += solver().height(lineSite, boundarySite);
        }

        else
        {
            hAvg += solver().height(boundarySite, lineSite);
        }
    }

    return hAvg/(double)depth();
}

void TrackLineAverage::initializeObserver(const Subjects &subject)
{
    (void) subject;
    m_averages.zeros();

    if (dim() == 0)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            m_averages(y) = bruteForceAverage(y);
        }
    }

    else
    {
        for (uint x = 0; x < solver().length(); ++x)
        {
            m_averages(x) = bruteForceAverage(x);
        }
    }

    onAverageChange();
}

void TrackLineAverage::notifyObserver(const Subjects &subject)
{
    (void) subject;

    const CurrentSurfaceChange &csc = solver().currentSurfaceChange();

    const uint &x0 = csc.x;
    const uint &y0 = csc.y;

    //if it is inside the area it is on some line.
    const bool initialInside = pointIsInsideArea(x0, y0);


    if (csc.type == ChangeTypes::Single)
    {
        if (initialInside)
        {
            if (dim() == 0)
            {
                m_averages(y0) += csc.value/(double)depth();
            }

            else
            {
                m_averages(x0) += csc.value/(double)depth();
            }

            onAverageChange();
        }
    }

    else
    {
        const int &x1 = csc.x1;
        const int &y1 = csc.y1;

        const bool finalInside = pointIsInsideArea(x1, y1);

        if (initialInside)
        {
            if (dim() == 0)
            {
                m_averages(y0) -= 1./depth();
            }

            else
            {
                m_averages(x0) -= 1./depth();
            }
        }

        if (finalInside)
        {
            if (dim() == 0)
            {
                m_averages(y1) += 1./depth();
            }

            else
            {
                m_averages(x1) += 1./depth();
            }
        }

        if (initialInside || finalInside)
        {
            onAverageChange();
        }
    }
}


double TrackAreaAverage::getBruteForceAverage() const
{
    int sum = 0;

    for (uint boundarySite = 0; boundarySite < solver().orthExtent(dim()); ++boundarySite)
    {
        for (uint lineSite = x0(); lineSite < x1(); ++lineSite)
        {
            if (dim() == 0)
            {
                sum += solver().height(lineSite, boundarySite);
            }

            else
            {
                sum += solver().height(boundarySite, lineSite);
            }
        }
    }

    return sum/(double)area();
}

void TrackAreaAverage::notifyObserver(const Subjects &subject)
{
    (void) subject;

    const CurrentSurfaceChange &csc = solver().currentSurfaceChange();
    const uint &x0 = csc.x;
    const uint &y0 = csc.y;
    const bool initialInside = pointIsInsideArea(x0, y0);

    if (csc.type == ChangeTypes::Single)
    {
        if (initialInside)
        {
            m_average += csc.value/(double)area();

            onAverageChange();
        }
    }

    else
    {
        const int &x1 = csc.x1;
        const int &y1 = csc.y1;

        const bool finalInside = pointIsInsideArea(x1, y1);

        if (initialInside && !finalInside)
        {
            m_average -= 1./area();

            onAverageChange();
        }

        else if (!initialInside && finalInside)
        {
            m_average += 1./area();

            onAverageChange();
        }

    }
}


double PartialBoundaryNeighbors::potential(const uint x, const uint y) const
{
    if (m_tracker.pointIsAtBoundary(x, y))
    {
        const int &h = solver().height(x, y);
        const double averageHeight = m_tracker.averageHeight(x, y);
        const int base = floor(averageHeight);

        double overlap;

        //There is no overlap between estimated height and current height
        if (h > ceil(averageHeight))
        {
            overlap = 0.0;
        }

        //There is 100% overlap between estimated and current height
        else if (h <= base)
        {
            overlap = 1.0;
        }

        else
        {
            overlap = averageHeight - base;
        }

        //We assume the boundary is blocked already, so we subtract one
        //to get rid of the neighbor interaction, and replace it by
        //the overlap, which is sort of a partial neighbor interaction.
        return overlap - 1.0;
    }

    return 0.0;
}


double NoBoundaryNeighbors::potential(const uint x, const uint y) const
{
    if (pointIsOnBoundary(x, y))
    {
        const int &h = solver().height(x, y);

        if (h > m_surfaceLevel)
        {
            return -1.;
        }
    }

    return 0.;
}
