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


double ConstantHeight::transformContinousCoordinate(const double xi, const double xj, const double xk) const
{
    return transformFunction<double>(xi, xj, xk);
}

bool ConstantHeight::isBlockedContinous(const double xi, const double xj, const double xk) const
{
    return blockedFunction<double>(xi, xj, xk);
}

bool ConstantHeight::isBlockedLattice(const int xi, const int xj, const int xk) const
{
    return blockedFunction<int>(xi, xj, xk);
}

int ConstantHeight::transformLatticeCoordinate(const int xi, const int xj, const int xk) const
{
    return transformFunction<int>(xi, xj, xk);
}

void ConstantHeight::closestImage(const double xi, const double xj, const double xk,
                                  const double xti, const double xtj, const double xtk,
                                  double &dxi, double &dxj, double &dxk) const
{
    noImage(xi, xj, xk, xti, xtj, xtk, dxi, dxj, dxk);
}
