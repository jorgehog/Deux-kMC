#pragma once

#include "solidonsolidevent.h"

class SolidOnSolidSystem;

class PressureWall : public SolidOnSolidEvent
{
public:

    PressureWall(SolidOnSolidSolver &solver,
                 const double E0,
                 const double sigma0,
                 const double r0);

    ~PressureWall();

    static double expSmallArg(double arg);

    void setupInitialConditions();

    void execute();

    void reset();

    const double &heightChange() const
    {
        return m_heightChange;
    }

    void recalculateRDLEnergy(const uint x, const uint y)
    {
        m_RDLEnergy(x, y) = EvaluateRDLEnergy(x, y);
    }

    double RDLEnergy(const uint x, const uint y) const
    {
        if (m_onsetTime != 0) BADAssClose(m_RDLEnergy(x, y), EvaluateRDLEnergy(x, y), 1E-5);

        return m_RDLEnergy(x, y);
    }

    double partialThetaRatio(const uint x, const uint y) const;

    double EvaluateRDLEnergy(const uint x, const uint y) const
    {
        const int &h = solver().height(x, y);
        return _RDLEnergyExpression(m_height - h);
    }

    const double &height() const
    {
        return m_height;
    }

    double RDLEnergySum() const;

    void findNewHeight();

    void updateRatesFor(DiffusionDeposition &reaction);

    void registerHeightChange(const uint x, const uint y, std::vector<DiffusionDeposition *> affectedReactions, const uint n);

    double bruteForceThetaRatio() const;

    const double &debyeLength() const
    {
        return m_lD;
    }

private:

    double m_height;
    double m_heightChange;
    double m_expFac;

    //    double m_thetaShift;

    //    long double m_thetaPrev;
    long double m_r0LogThetaPrev;

    const double m_lD;
    const double m_s0;
    const double m_E0;

    mat m_RDLEnergy;
    //    vector<vector<long double>> m_thetaPartialSumsMat;
    mat m_ratioPartialSums;

    void setupTheta();

    void recalculateAllRDLEnergies();

    double _RDLEnergyExpression(const double heightDifference) const
    {
        return -m_s0*std::exp(-heightDifference/m_lD);
    }

};
