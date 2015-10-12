#include "constantheight.h"

using namespace kMC;

ConstantHeight::ConstantHeight(const double height,
                               const double location,
                               const Boundary::orientations orientation) :
    Boundary(orientation),
    m_height(height),
    m_location(location)
{

}


ConstantHeight::~ConstantHeight()
{

}


double ConstantHeight::transformCoordinate(const double xi, const double xj, const double xk) const
{
    (void) xj;
    (void) xk;

    return xi;
}

bool ConstantHeight::isBlocked(const double xi, const double xj, const double xk) const
{
    (void) xj;

    if (orientation() == orientations::FIRST)
    {
        return (xi < m_location) && (xk <= m_height);
    }

    else
    {
        return (xi > m_location) && (xk <= m_height);
    }
}

void ConstantHeight::closestImage(const double xi, const double xj, const double xk,
                                  const double xti, const double xtj, const double xtk,
                                  double &dxi, double &dxj, double &dxk) const
{
    noImage(xi, xj, xk, xti, xtj, xtk, dxi, dxj, dxk);
}
