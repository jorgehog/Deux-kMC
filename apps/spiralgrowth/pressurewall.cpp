#include "pressurewall.h"

#include "solidonsolidreaction.h"


PressureWall::PressureWall(SolidOnSolidSolver &solver, const double E0,
                           const double sigma0,
                           const double r0):
    SolidOnSolidEvent(solver, "PressureWall", "h0", true, true),
    m_r0(r0),
    m_s0(sigma0),
    m_E0(E0),
    m_localPressure(solver.length(), solver.width(), fill::zeros),
    m_thetaPartialSumsMat(solver.length(), solver.width())
{
    BADAss(E0, >, 0, "E0 should be positive.");

    solver.setPressureWallEvent(*this);
 }

PressureWall::~PressureWall()
{

}

double PressureWall::expSmallArg(double arg)
{
    if (arg > 0.1 || arg < -0.1)
    {
        return std::exp(arg);
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

void PressureWall::setupInitialConditions()
{
    setupTheta();

    m_height = m_r0*std::log(solver().area()*m_thetaPrev/(m_E0/m_s0));

    recalculateAllPressures();

    BADAssClose(pressureEnergySum(), -m_E0, 1E-5);

}

void PressureWall::execute()
{
    double value = m_height - dependency("AverageHeight")->value();

    setValue(value);


    if ((cycle() + 1)%100 == 0)
    {
        recalculateAllPressures();
    }
}

void PressureWall::reset()
{

    //Check that everything is updated correctly from previous runs.
#ifndef NDEBUG
    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            BADAssClose(localPressureEvaluate(x, y), localPressure(x, y), 1E-3);
        }
    }
#endif

    BADAssClose(pressureEnergySum(), -m_E0, 1E-5);

}

double PressureWall::partialTheta(const uint x, const uint y) const
{
    return std::exp((solver().height(x, y))/m_r0);
}

double PressureWall::pressureEnergySum() const
{
    double s = 0;
    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            BADAssClose(m_localPressure(x, y), localPressureEvaluate(x, y), 1E-4);

            s += m_localPressure(x, y);
        }
    }

    return s;
}


void PressureWall::findNewHeight()
{

    double theta = accu(m_thetaPartialSumsMat/solver().area());

//    BADAssClose(theta, bruteForceTheta(), 1E-5);

//    theta = bruteForceTheta();

    m_heightChange = m_r0*std::log(theta/m_thetaPrev);

    m_thetaPrev = theta;

    m_height += m_heightChange;

    m_expFac = expSmallArg(-m_heightChange/m_r0);

}

void PressureWall::updateRatesFor(DiffusionDeposition &reaction)
{
    const uint &x = reaction.x();
    const uint &y = reaction.y();

    double rateChange = expSmallArg(-solver().alpha()*m_localPressure(x, y)*(m_expFac - 1));

    //For every affected particle we update only those who include the pressure term.
    //Vector is set up in initialize based on virtual reaction function isPressureAffected().

    double prevDiffRate = reaction.diffusionRate();
    reaction.setDiffusionRate(prevDiffRate*rateChange);
    reaction.changeRate(reaction.diffusionRate() + reaction.depositionRate());

    m_localPressure(x, y) *= m_expFac;

    BADAssClose(localPressureEvaluate(x, y), m_localPressure(x, y), 1E-3, "incorrect pressure update", [&] ()
    {
        BADAssSimpleDump(cycle(), x, y, localPressure(x, y), m_expFac, m_heightChange);
    });

    BADAssClose(reaction.diffusionRate(), reaction.calculateDiffusionRate(), 1E-5, "incorrect rate update", [&] ()
    {
        BADAssSimpleDump(cycle(), x, y, localPressure(x, y), m_expFac, m_heightChange);
    });

}

void PressureWall::registerHeightChange(const uint x, const uint y)
{
    m_thetaPartialSumsMat(x, y) = partialTheta(x, y);
}

double PressureWall::bruteForceTheta() const
{
    double theta = 0;

    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            theta += partialTheta(x, y);
        }
    }

    return theta/solver().area();
}

void PressureWall::setupTheta()
{
    m_thetaPrev = 0;

    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            m_thetaPartialSumsMat(x, y) = partialTheta(x, y);
            m_thetaPrev += m_thetaPartialSumsMat(x, y);
        }
    }

    m_thetaPrev /= solver().area();

}

void PressureWall::recalculateAllPressures()
{
    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            m_localPressure(x, y) = localPressureEvaluate(x, y);
        }
    }
}

