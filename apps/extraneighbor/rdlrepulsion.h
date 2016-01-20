#pragma once

#include <SOSkMC.h>


class RDLRepulsion : public LocalCachedPotential
{
public:
    RDLRepulsion(SOSSolver &solver, const double s0, const double ld, const double Pl);

private:
    const double m_s0;
    const double m_ld;
    const double m_Pl;

    // Observer interface
public:
    void notifyObserver(const Subjects &subject);

    // LocalCachedPotential interface
public:
    double potentialFunction(const uint x, const uint y) const;
};

