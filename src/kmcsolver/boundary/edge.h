#pragma once

#include "finitesize.h"

namespace kMC
{

class Edge : public FiniteSize
{
public:

    Edge(const int location,
         const orientations orientation);

    bool isBlocked(const double xi) const;

private:

    const int m_location;

};

}
