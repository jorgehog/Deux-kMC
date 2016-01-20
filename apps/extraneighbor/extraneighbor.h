#pragma once

#include <SOSkMC.h>

class ExtraNeighbor : public LocalCachedPotential
{
public:
    ExtraNeighbor(SOSSolver &solver);

    double energyFunction(const double dh) const;

    static constexpr double m_shift = pow(2, -6.);
    static constexpr double m_scaling = 1./(1 - m_shift);

    // Observer interface
public:
    void notifyObserver(const Subjects &subject);

    // LocalCachedPotential interface
public:
    double potentialFunction(const uint x, const uint y) const;
};

