#pragma once

#include "boundary.h"

namespace kMC
{

class Edge : public Boundary
{
public:

    Edge(const uint span);

    int transformCoordinate(const int xi) const;

    bool isBlocked(const int xi) const;

private:

    const uint m_span;

};

}
