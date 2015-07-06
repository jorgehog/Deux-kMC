#ifndef RDLSURFACE_H
#define RDLSURFACE_H

#include "confiningsurface.h"

class RDLSurface : public virtual ConfiningSurface
{
public:
    RDLSurface(SolidOnSolidSolver &solver,
               const double E0,
               const double s0,
               const double lD);

    virtual ~RDLSurface();

    static double expSmallArg(double arg);

    double partialThetaRatio(const uint x, const uint y) const;

    double RDLEnergySum() const;

    void findNewHeight();

    void updateRatesFor(DiffusionDeposition &reaction);

    double bruteForceThetaRatio() const;

    void recalculateAllRDLEnergies();

    const double &debyeLength() const
    {
        return m_lD;
    }

    double evaluateRDLEnergy(const uint x, const uint y) const
    {
        return _RDLEnergyExpression(height() - solver().height(x, y));
    }

    double RDLEnergy(const uint x, const uint y) const
    {
        if (m_onsetTime != 0) BADAssClose(m_RDLEnergy(x, y), evaluateRDLEnergy(x, y), 1E-5);

        return m_RDLEnergy(x, y);
    }

    double _RDLEnergyExpression(const double heightDifference) const
    {
        return -m_s0*std::exp(-heightDifference/m_lD);
    }

    const double &heightChange() const
    {
        return m_heightChange;
    }

    void recalculateRDLEnergy(const uint x, const uint y)
    {
        m_RDLEnergy(x, y) = evaluateRDLEnergy(x, y);
    }


private:

    double m_heightChange;
    double m_expFac;

    long double m_r0LogThetaPrev;

    const double m_lD;
    const double m_s0;
    const double m_E0;

    mat m_RDLEnergy;
    mat m_ratioPartialSums;

    void setupTheta();


    // Event interface
public:
    void execute();

    void initialize();

    void reset();

    // ConfiningSurface interface
public:
    void setupInitialConditions();

    void registerHeightChange(const uint x, const uint y, std::vector<DiffusionDeposition *> affectedReactions, const uint n);

    double confinementEnergy(const uint x, const uint y)
    {
        return RDLEnergy(x, y);
    }
};

#endif // RDLSURFACE_H
