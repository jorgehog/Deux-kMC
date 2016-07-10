#pragma once

#include <SOSkMC.h>

class ExtraNeighbor : public LocalCachedPotential
{
public:
    ExtraNeighbor(SOSSolver &solver, const double relBondEnergy = 1.0);

    double energyFunction(const double dh) const;

    const double &relBondEnergy() const
    {
        return m_relBondEnergy;
    }

private:
    const double m_relBondEnergy;

    // Observer interface
public:
    void notifyObserver(const Subjects &subject);

    // LocalCachedPotential interface
public:
    double potentialFunction(const uint x, const uint y) const;
};

