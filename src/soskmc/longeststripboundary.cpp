#include "longeststripboundary.h"

LongestStripBoundary::LongestStripBoundary(SOSSolver &solver,
                                           const uint dim,
                                           orientations orientation) :
    SOSBoundary(solver, orientation),
    m_dim(dim)
{
    if (dim == 0)
    {
        m_yspan = solver.width();

        if (orientation == Boundary::orientations::FIRST)
        {
            m_location = 0;
        }

        else
        {
            m_location = solver.length() - 1;
        }
    }

    else
    {
        m_yspan = solver.length();

        if (orientation == Boundary::orientations::FIRST)
        {
            m_location = 0;
        }

        else
        {
            m_location = solver.width() - 1;
        }
    }

    solver.registerObserver(this);
}


bool LongestStripBoundary::isConnectedToOriginCalc(const uint &y, const int &h) const
{
    if (h < m_currentHeight)
    {
        return true;
    }

    else if (h > m_currentHeight)
    {
        return false;
    }

    if (m_dim == 0)
    {
        if (orientation() == orientations::FIRST)
        {
            for (uint ys = 0; ys < y; ++ys)
            {
                if (solver().height(0, ys) < h)
                {
                    return false;
                }
            }
        }

        else
        {
            for (uint ys = 0; ys < y; ++ys)
            {
                if (solver().height(solver().length() - 1, ys) < h)
                {
                    return false;
                }
            }
        }
    }

    else
    {
        if (orientation() == orientations::FIRST)
        {
            for (uint ys = 0; ys < y; ++ys)
            {
                if (solver().height(ys, 0) < h)
                {
                    return false;
                }
            }
        }

        else
        {
            for (uint ys = 0; ys < y; ++ys)
            {
                if (solver().height(ys, solver().width() - 1) < h)
                {
                    return false;
                }
            }
        }
    }

    return true;
}

bool LongestStripBoundary::isConnectedToOrigin(const uint y, const int h) const
{
    return (h < m_currentHeight) || ((h == m_currentHeight) && (y <= m_ymax || y >= m_ymax2));
}

void LongestStripBoundary::scanForYMax()
{
    if (m_dim == 0)
    {
        while (solver().height(0, m_ymax+1) >= m_currentHeight)
        {
            m_ymax++;

            if (m_ymax == solver().width() - 1)
            {
                m_ymax = 0;
                m_currentHeight++;
                return scanForYMax();
            }
        }
    }

    else
    {
        while (solver().height(m_ymax+1, 0) >= m_currentHeight)
        {
            m_ymax++;

            if (m_ymax == solver().length() - 1)
            {
                m_ymax = 0;
                m_currentHeight++;
                return scanForYMax();
            }
        }
    }
}

void LongestStripBoundary::scanForYMax2()
{
    if (m_dim == 0)
    {
        while (solver().height(0, m_ymax2-1) >= m_currentHeight)
        {
            m_ymax2--;

            if (m_ymax2 == 0)
            {
                m_ymax2 = solver().width() - 1;
                m_currentHeight++;
                return scanForYMax2();
            }
        }
    }

    else
    {
        while (solver().height(m_ymax2-1, 0) >= m_currentHeight)
        {
            m_ymax2--;

            if (m_ymax2 == 0)
            {
                m_ymax2 = solver().length() - 1;
                m_currentHeight++;
                return scanForYMax2();
            }
        }
    }
}

double LongestStripBoundary::transformContinousCoordinate(const double xi, const double xj, const double xk) const
{
    (void) xj;
    (void) xk;

    return xi;
}

int LongestStripBoundary::transformLatticeCoordinate(const int xi, const int xj, const int xk) const
{
    (void) xj;
    (void) xk;

    return xi;
}

bool LongestStripBoundary::isBlockedContinous(const double xi, const double xj, const double xk) const
{
    (void) xi;
    (void) xj;
    (void) xk;

    return false;
}

bool LongestStripBoundary::isBlockedLattice(const int xi, const int xj, const int xk) const
{
    if (orientation() == orientations::FIRST)
    {
        if (xi < m_location)
        {
            return isConnectedToOrigin(xj, xk);
        }
    }

    else
    {
        if (xi >= m_location)
        {
            return isConnectedToOrigin(xj, xk);
        }
    }

    return false;
}

void LongestStripBoundary::closestImage(const double xi, const double xj, const double xk, const double xti, const double xtj, const double xtk, double &dxi, double &dxj, double &dxk) const
{
    noImage(xi, xj, xk, xti, xtj, xtk, dxi, dxj, dxk);
}

void LongestStripBoundary::initializeObserver(const Subjects &subject)
{
    (void) subject;

    m_currentHeight = solver().heights().min();

    m_ymax = 0;
    m_ymax2 = m_dim == 0 ? solver().width() : solver().length();
    scanForYMax();
    scanForYMax2();

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

void LongestStripBoundary::notifyObserver(const Subjects &subject)
{
    (void) subject;

    const uint ymaxPrev = m_ymax;
    const uint ymax2Prev = m_ymax2;
    scanForYMax();
    scanForYMax2();

    if (ymaxPrev == m_ymax && ymax2Prev == m_ymax2)
    {
        return;
    }

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
