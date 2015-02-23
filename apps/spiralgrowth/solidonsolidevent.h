#pragma once

#include <kMC>
#include <utils.h>
#include "solidonsolidsolver.h"

using namespace ignis;

class SolidOnSolidEvent : public ignis::LatticeEvent
{
public:
    using LatticeEvent::LatticeEvent;
    virtual ~SolidOnSolidEvent() {}

    const SolidOnSolidSolver *solver() const
    {
        return dependency<SolidOnSolidSolver>("KMCSolver");
    }

};
