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


double kMC::ConstantHeight::transformCoordinate(const double xi, const double xj, const double xk) const
{
    (void) xj;
    (void) xk;

    return xi;
}

bool kMC::ConstantHeight::isBlocked(const double xi, const double xj, const double xk) const
{
    (void) xj;

    if (m_orientation == orientations::FIRST)
    {
        return (xi < m_location) && (xk <= m_height);
    }

    else
    {
        return (xi > m_location) && (xk <= m_height);
    }
}

std::vector<double> kMC::ConstantHeight::imagesOf(const double xi, const double xj, const double xk) const
{
    (void) xi;
    (void) xj;
    (void) xk;

    return {};
}
