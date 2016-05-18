#pragma once

#include <SOSkMC.h>

class ExtraNeighbor : public LocalCachedPotential
{
public:
    ExtraNeighbor(SOSSolver &solver);

    double energyFunction(const double dh) const;

    // Observer interface
public:
    void notifyObserver(const Subjects &subject);

    // LocalCachedPotential interface
public:
    double potentialFunction(const uint x, const uint y) const;
};

