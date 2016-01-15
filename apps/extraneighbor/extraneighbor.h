#pragma once

#include <SOSkMC.h>

class ExtraNeighbor : public LocalPotential, public Observer<Subjects>
{
public:
    ExtraNeighbor(SOSSolver &solver);

    double energyFunction(const double dh) const;
    double energyFunction(const uint x, const uint y) const;

private:
    mat m_potentialValues;

    // Observer interface
public:
    void initializeObserver(const Subjects &subject);
    void notifyObserver(const Subjects &subject);

    // LocalPotential interface
public:
    double energy(const uint x, const uint y) const;
};

