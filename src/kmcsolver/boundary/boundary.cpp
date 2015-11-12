#include "boundary.h"

using namespace kMC;

Boundary::Boundary(const Boundary::orientations orientation) :
    m_orientation(orientation)
{

}



void Boundary::noImage(const double xi, const double xj, const double xk,
                       const double xti, const double xtj, const double xtk,
                       double &dxi, double &dxj, double &dxk) const
{
    (void) xj;
    (void) xk;
    (void) xtj;
    (void) xtk;

    dxi = xti - xi;
    dxj = std::numeric_limits<double>::max();
    dxk = std::numeric_limits<double>::max();
}
