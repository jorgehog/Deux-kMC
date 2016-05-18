#pragma once

#include <SOSkMC.h>

class ExtraNeighbor : public LocalCachedPotential
{
public:
    ExtraNeighbor(SOSSolver &solver, const double energyShift);

    double energyFunction(const double dh) const;

private:

    const double m_energyShift;

    // Observer interface
public:
    void notifyObserver(const Subjects &subject);

    // LocalCachedPotential interface
public:
    double potentialFunction(const uint x, const uint y) const;
};

