#pragma once

#include "confiningsurface.h"

class SurfaceReaction;
class RDLPotential;

class RDLSurface : public virtual ConfiningSurface
{
public:
    RDLSurface(SOSSolver &solver, RDLPotential &potential, const double Pl);

    virtual ~RDLSurface();

    double partialThetaRatio(const uint x, const uint y) const;

    double RDLEnergySum() const;

    void findNewHeight();

    double bruteForceThetaRatio() const;

private:

    const RDLPotential &m_potential;

    const double m_Pl;

    long double m_ldLogThetaPrev;

    mat m_ratioPartialSums;

    void setupTheta();

    double calculateKZrel(const double x, const double y, const double z, double &K, double &zRel) const;

    double calculateKZrel(const double x, const double y, const double z, double &K) const
    {
        double zRel;
        return calculateKZrel(x, y, z, K, zRel);
    }

    double calculateKZrel(const double x, const double y, const double z) const
    {
        double K, zRel;
        return calculateKZrel(x, y, z, K, zRel);
    }

    // Event interface
public:
    void execute();
    void initialize();
    void reset();

    // ConfiningSurface interface
public:

    bool hasSurface() const
    {
        return true;
    }

    bool acceptDiffusionMove(const double x0, const double y0, const double z0, const double x1, const double y1, const double z1) const;

    double diffusionDrift(const double x, const double y, const double z) const;


    // kMC::Observer interface
public:

    void initializeObserver(const Subjects &subject);

    void notifyObserver(const Subjects &subject);
};
