#include "averageheightlineboundary.h"

AverageHeightLineBoundary::AverageHeightLineBoundary(SOSSolver &solver,
                                                     const Boundary::orientations orientation,
                                                     const uint dim,
                                                     const uint depth) :
    SOSBoundary(solver, orientation),
    m_dim(dim),
    m_depth(depth)
{
    if (orientation == Boundary::orientations::FIRST)
    {
        m_x0 = 0;
        m_x1 = depth;
    }

    else
    {
        if (dim == 0)
        {
            m_x0 = solver.length() - depth;
            m_x1 = solver.length();
        }

        else
        {
            m_x0 = solver.width() - depth;
            m_x1 = solver.width();
        }
    }
}

double AverageHeightLineBoundary::transformContinousCoordinate(const double xi, const double xj, const double xk) const
{
    (void) xj;
    (void) xk;
    return xi;
}

int AverageHeightLineBoundary::transformLatticeCoordinate(const int xi, const int xj, const int xk) const
{
    (void) xj;
    (void) xk;
    return xi;
}

bool AverageHeightLineBoundary::isBlockedContinous(const double xi, const double xj, const double xk) const
{
    (void) xi;
    (void) xj;
    (void) xk;

    return false;
}

bool AverageHeightLineBoundary::isBlockedLattice(const int xi, const int xj, const int xk) const
{
    if (orientation() == Boundary::orientations::FIRST)
    {
        if (xi >= 0)
        {
            return false;
        }
    }

    else
    {
        if (m_dim == 0)
        {
            if (xi < int(solver().length()))
            {
                return false;
            }
        }

        else
        {
            if (xi < int(solver().width()))
            {
                return false;
            }
        }
    }

    double hAvg = 0;

    if (m_dim == 0)
    {
        for (uint x = m_x0; x < m_x1; ++x)
        {
            hAvg += solver().height(x, xj);
        }

        hAvg /= m_depth;

        if (orientation() == Boundary::orientations::FIRST)
        {
            if (solver().height(0, xj) < round(hAvg))
            {
                return false;
            }

            else
            {
                return true;
            }
        }

        else
        {
            if (solver().height(solver().length() - 1, xj) < round(hAvg))
            {
                return false;
            }

            else
            {
                return true;
            }
        }
    }

    else
    {
        for (uint y = m_x0; y < m_x1; ++y)
        {
            hAvg += solver().height(xj, y);
        }

        hAvg /= m_depth;

        if (orientation() == Boundary::orientations::FIRST)
        {
            if (xk <= round(hAvg))
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
            if (xk <= round(hAvg))
            {
                return true;
            }

            else
            {
                return false;
            }
        }

    }

}

void AverageHeightLineBoundary::closestImage(const double xi, const double xj, const double xk, const double xti, const double xtj, const double xtk, double &dxi, double &dxj, double &dxk) const
{
    return noImage(xi, xj, xk, xti, xtj, xtk, dxi, dxj, dxk);
}
