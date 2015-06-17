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

    void recalculateLocalPressure(const uint x, const uint y)
    {
        m_localPressure(x, y) = localPressureEvaluate(x, y);
    }

    double localPressure(const uint x, const uint y) const
    {
        if (m_onsetTime != 0) BADAssClose(m_localPressure(x, y), localPressureEvaluate(x, y), 1E-5);

        return m_localPressure(x, y);
    }

    double partialRatio(const uint x, const uint y) const;

    double localPressureEvaluate(const uint x, const uint y) const
    {
        const int &h = solver().height(x, y);
        return _pressureExpression(m_height - h);
    }

    const double &height() const
    {
        return m_height;
    }

    double pressureEnergySum() const;

    void findNewHeight();

    void updateRatesFor(DiffusionDeposition &reaction);

    void registerHeightChange(const uint x, const uint y, std::vector<DiffusionDeposition *> affectedReactions, const uint n);

    double bruteForceRatio() const;

    const double &debyeLength() const
    {
        return m_r0;
    }

private:

    double m_height;
    double m_heightChange;
    double m_expFac;

    //    double m_thetaShift;

    //    long double m_thetaPrev;
    long double m_r0LogThetaPrev;

    const double m_r0;
    const double m_s0;
    const double m_E0;

    mat m_localPressure;
    //    vector<vector<long double>> m_thetaPartialSumsMat;
    mat m_ratioPartialSums;

    void setupTheta();

    void recalculateAllPressures();

    double _pressureExpression(const double heightDifference) const
    {
        return -m_s0*std::exp(-heightDifference/m_r0);
    }

};
