#include "miscboundaries.h"

AverageHeightLineBoundary::AverageHeightLineBoundary(SOSSolver &solver,
                                                     const Boundary::orientations orientation,
                                                     const uint dim,
                                                     const uint depth) :
    SOSBoundary(solver, orientation),
    m_dim(dim),
    m_depth(depth)
{
    solver.registerObserver(this);

    if (orientation == Boundary::orientations::FIRST)
    {
        m_x0 = 0;
        m_x1 = depth;
        m_location = 0;
    }

    else
    {
        if (dim == 0)
        {
            m_x0 = solver.length() - depth;
            m_x1 = solver.length();
            m_location = solver.length();
        }

        else
        {
            m_x0 = solver.width() - depth;
            m_x1 = solver.width();
            m_location = solver.width();
        }
    }
}

int AverageHeightLineBoundary::transformLatticeCoordinate(const int xi, const int xj, const int xk) const
{
    (void) xj;
    (void) xk;
    return xi;
}

bool AverageHeightLineBoundary::isBlockedLattice(const int xi, const int xj, const int xk) const
{
    if (orientation() == First)
    {
        if (xi >= 0)
        {
            return false;
        }
    }

    else
    {
        if (xi < int(m_location))
        {
            return false;
        }
    }

    double hAvg = 0;

    for (uint coor = m_x0; coor < m_x1; ++coor)
    {
        if (m_dim == 0)
        {
            hAvg += solver().height(coor, xj);
        }

        else
        {
            hAvg += solver().height(xj, coor);
        }
    }

    hAvg /= m_depth;

    const int hBoundary = int(std::round(hAvg));

    if (orientation() == First)
    {
        if (xk <= hBoundary)
        {
            return true;
        }

        else
        {
            return false;
        }
    }

    else
    {
        if (xk <= hBoundary)
        {
            return true;
        }

        else
        {
            return false;
        }
    }

}

void AverageHeightLineBoundary::closestImage(const double xi, const double xj, const double xk, const double xti, const double xtj, const double xtk, double &dxi, double &dxj, double &dxk) const
{
    return noImage(xi, xj, xk, xti, xtj, xtk, dxi, dxj, dxk);
}


void AverageHeightLineBoundary::initializeObserver(const Subjects &subject)
{
    (void) subject;
}

void AverageHeightLineBoundary::notifyObserver(const Subjects &subject)
{
    (void) subject;

    uint loc;
    if (m_dim == 0)
    {
        if (orientation() == First)
        {
            loc = 0;
        }

        else
        {
            loc = solver().length() - 1;
        }

        for (uint y = 0; y < solver().width(); ++y)
        {
            solver().registerChangedSite(loc, y);
        }
    }

    else
    {
        if (orientation() == First)
        {
            loc = 0;
        }

        else
        {
            loc = solver().width() - 1;
        }

        for (uint x = 0; x < solver().length(); ++x)
        {
            solver().registerChangedSite(x, loc);
        }
    }
}


ReflConstantHybrid::ReflConstantHybrid(SOSSolver &solver, const int h, const Boundary::orientations orientation, const uint dim) :
    SOSBoundary(solver, orientation),
    m_h(h)
{
    if (orientation == First)
    {
        m_location = 0;
    }

    else
    {
        if (dim == 0)
        {
            m_location = solver.length();
        }

        else
        {
            m_location = solver.width();
        }
    }
}
