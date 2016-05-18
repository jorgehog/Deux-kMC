#pragma once

#include "../../localcachedpotential.h"

class RDLSurface;

class RDLPotential : public LocalCachedPotential
{
public:
    RDLPotential(SOSSolver &solver, const double s0, const double lD);

    const double &expFac() const
    {
        return m_expFac;
    }

    const double &lD() const
    {
        return m_lD;
    }

    const double &s0() const
    {
        return m_s0;
    }

    double rdlEnergy(const double dh) const;

    double rdlEnergyDeriv(const double dh) const;

    static double expSmallArg(double arg);

private:

    double m_expFac;

    const double m_s0;
    const double m_lD;

    // Observer interface
public:
    void notifyObserver(const Subjects &subject);

    // LocalCachedPotential interface
public:
    double potentialFunction(const uint x, const uint y) const;
};
