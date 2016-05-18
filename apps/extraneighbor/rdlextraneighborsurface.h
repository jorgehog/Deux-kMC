#pragma once

#include <SOSkMC.h>
#include "extraneighbor.h"

class RDLExtraNeighborSurface : public ConfiningSurface
{
public:
    RDLExtraNeighborSurface(SOSSolver &solver,
                            RDLPotential &rdlPotential,
                            ExtraNeighbor &extraNeighborPotential, const double Pl);

    double getRdlEquilibrium() const;

    void dumpProfile() const;

    double totalForce(const double hl) const;

    double totalForceDeriv(const double hl) const;

    double totalWdVAttraction(const double hl) const;

private:

    const RDLPotential &m_rdlPotential;
    const ExtraNeighbor &m_extraNeighborPotential;

    const double m_Pl;
    const double m_nFactor;

    void findNewHeight();

    double getHeightBisection() const;
    double bisect(const double min, const double max, double fmin,
                  const double eps=1E-12,
                  const uint nMax = std::numeric_limits<uint>::max()) const;

    // Observer interface
public:
    void initializeObserver(const Subjects &subject);
    void notifyObserver(const Subjects &subject);

    // Event interface
public:
    void execute();

    // ConfiningSurface interface
public:
    bool hasSurface() const;
    bool acceptDiffusionMove(const double x0, const double y0, const double z0, const double x1, const double y1, const double z1) const;
    double diffusionDrift(const double x, const double y, const double z) const;
};
