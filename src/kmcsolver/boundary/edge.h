#pragma once

#include "finitesize.h"

namespace kMC
{

class Edge : public FiniteSize
{
public:

    Edge(const double location,
         const orientations orientation);

    bool isBlocked(const double xi) const;

private:

    const double m_location;

};

}
