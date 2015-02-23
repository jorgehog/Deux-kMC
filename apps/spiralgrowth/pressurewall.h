#pragma once

#include "solidonsolidevent.h"

namespace kMC
{

class SolidOnSolidSystem;

class PressureWall : public SolidOnSolidEvent
{
public:

    PressureWall(SolidOnSolidSolver &mutexSolver, const double E0,
               const double sigma0,
               const double r0);

    const string numericDescription() const
    {
        stringstream s;
        s << "E0_" << m_E0
          << "_s0_" << m_s0
          << "_r0_" << m_r0;

        return s.str();

    }

    ~PressureWall();

    static double expSmallArg(double arg)
    {
        if (arg > 0.1 || arg < -0.1)
        {
            return exp(arg);
        }

        BADAssClose(arg, 0, 0.1, "Argument is not small.", [&arg] ()
        {
            BADAssSimpleDump(arg);
        });

        double arg2 = arg*arg;
        double arg4 = arg2*arg2;
        double approx = 1.0 + arg*(1 + 1.0/6*arg2 + 1.0/120*arg4) + 0.5*(arg2 + 1.0/12*arg4);

        BADAssClose(exp(arg), approx, 1E-5,
                    "Exponential approximation failed.", [&] ()
        {
            BADAssSimpleDump(arg, exp(arg), approx);
        });

        return approx;
    }

    void initialize();

    void execute()
    {

    }

    void reset();

    const double &heightChange() const
    {
        return m_heightChange;
    }

    double localPressure(const uint x, const uint y) const
    {
        if (m_onsetTime != 0) BADAssClose(m_localPressure(x, y), localPressureEvaluate(x, y), 1E-5);

        return m_localPressure(x, y);
    }

    double localPressureEvaluate(const uint x, const uint y) const
    {
        return _pressureExpression(m_height - solver()->height(x, y));
    }

    const double &height() const
    {
        return m_height;
    }

    double pressureEnergySum() const
    {
        double s = 0;
        for (uint x = 0; x < solver()->length(); ++x)
        {
            for (uint y = 0; y < solver()->width(); ++y)
            {
                BADAssClose(m_localPressure(x, y), localPressureEvaluate(x, y), 1E-4);

                s += m_localPressure(x, y);
            }
        }

        return s;
    }

private:

    SolidOnSolidSolver &m_mutexSolver;

    double m_height;
    double m_heightChange;
    double m_thetaPrev;

    uint m_changed;

    const double m_r0;
    const double m_s0;
    const double m_E0;

    mat m_localPressure;

    SolidOnSolidSolver *solver() const
    {
        return &m_mutexSolver;
    }

    void _rescaleHeight();

    void _updatePressureRates();

    void recalculateAllPressures();

    double _pressureExpression(const double heightDifference) const
    {
        return -m_s0*std::exp(-heightDifference/m_r0);
    }

};


}
