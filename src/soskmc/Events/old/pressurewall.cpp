#include "pressurewall.h"

#include "solidonsolidreaction.h"


PressureWall::PressureWall(SolidOnSolidSolver &solver,
                           const double E0,
                           const double sigma0,
                           const double r0):
    SolidOnSolidEvent(solver, "PressureWall", "h0", true, true),
    m_lD(r0),
    m_s0(sigma0),
    m_E0(E0),
    m_RDLEnergy(solver.length(), solver.width(), fill::zeros),
    m_ratioPartialSums(solver.length(), solver.width())
{
    BADAss(E0, >, 0, "E0 should be positive.");

    solver.setConfiningSurfaceEvent(*this);
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

    m_height = m_r0LogThetaPrev + m_lD*std::log(m_s0/m_E0);

    recalculateAllRDLEnergies();

    BADAssClose(RDLEnergySum(), -m_E0, 1E-5);
}

void PressureWall::execute()
{
    double value = m_height - dependency("AverageHeight")->value();

    setValue(value);

    if ((cycle() + 1)%100 == 0)
    {
        recalculateAllRDLEnergies();
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
            BADAssClose(localPressureEvaluate(x, y), localPressure(x, y), 1E-5);
        }
    }
#endif

    BADAssClose(RDLEnergySum(), -m_E0, 1E-3);
}

double PressureWall::partialThetaRatio(const uint x, const uint y) const
{
    return std::exp((solver().height(x, y) - m_r0LogThetaPrev)/m_lD);
}

double PressureWall::RDLEnergySum() const
{
    double s = 0;
    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            BADAssClose(m_RDLEnergy(x, y), EvaluateRDLEnergy(x, y), 1E-5);

            s += m_RDLEnergy(x, y);
        }
    }

    return s;
}


void PressureWall::findNewHeight()
{
    double ratio = 0;

    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            ratio += m_ratioPartialSums(x, y);
            BADAssClose(m_ratioPartialSums(x, y), partialThetaRatio(x, y), 1E-3, "theta fail", [&]
            {
                BADAssSimpleDump(m_ratioPartialSums(x, y), partialThetaRatio(x, y), m_r0LogThetaPrev);
            });
        }
    }

    m_heightChange = m_lD*std::log(ratio);

    m_r0LogThetaPrev += m_heightChange;

    m_height += m_heightChange;

    m_expFac = expSmallArg(-m_heightChange/m_lD);

    m_ratioPartialSums*=m_expFac;
}

void PressureWall::updateRatesFor(DiffusionDeposition &reaction)
{
    const uint x = reaction.x();
    const uint y = reaction.y();

    double rateChange = expSmallArg(-solver().alpha()*m_RDLEnergy(x, y)*(m_expFac - 1));

    //For every affected particle we update only those who include the pressure term.
    //Vector is set up in initialize based on virtual reaction function isPressureAffected().

    double prevDiffRate = reaction.dissolutionRate();
    reaction.setDiffusionRate(prevDiffRate*rateChange);

    m_RDLEnergy(x, y) *= m_expFac;

    BADAssClose(EvaluateRDLEnergy(x, y), m_RDLEnergy(x, y), 1E-5, "incorrect pressure update", [&] ()
    {
        BADAssSimpleDump(cycle(), x, y, rateChange, localPressure(x, y), m_expFac, m_heightChange);
    });

    BADAssClose(reaction.dissolutionRate(), reaction.calculateDissolutionRate(), 1E-5, "incorrect rate update", [&] ()
    {
        double lp = localPressure(x, y);
        int h = solver().height(x, y);
        BADAssSimpleDump(cycle(), x, y, rateChange, prevDiffRate, lp, m_expFac, m_heightChange, h, m_height);
    });

}

void PressureWall::registerHeightChange(const uint x,
                                        const uint y,
                                        std::vector<DiffusionDeposition*> affectedReactions,
                                        const uint n)
{
    if (!hasStarted())
    {
        return;
    }

    auto start = affectedReactions.begin();
    auto end = start + n;

    m_ratioPartialSums(x, y) = partialThetaRatio(x, y);

    findNewHeight();

    for (uint i = 0; i < n; ++i)
    {
        DiffusionDeposition *reaction = affectedReactions.at(i);
        recalculateRDLEnergy(reaction->x(), reaction->y());
    }

    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            DiffusionDeposition &reaction = solver().reaction(x, y);

            if (std::find(start, end, &reaction) == end)
            {
                updateRatesFor(reaction);
            }
        }
    }
}


double PressureWall::bruteForceThetaRatio() const
{
    double ratio = 0;

    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            ratio += std::exp(((long double)solver().height(x, y) - m_r0LogThetaPrev)/m_lD);
        }
    }

    return ratio;
}

void PressureWall::setupTheta()
{
    double thetaShift = dependency("AverageHeight")->value();
    long double thetaPrev = 0;

    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            thetaPrev += std::exp((solver().height(x, y) - thetaShift)/m_lD);
        }
    }

    m_r0LogThetaPrev = m_lD*log(thetaPrev) + thetaShift;


    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            m_ratioPartialSums(x, y) = partialThetaRatio(x, y);
        }
    }

}

void PressureWall::recalculateAllRDLEnergies()
{
    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            m_RDLEnergy(x, y) = EvaluateRDLEnergy(x, y);
        }
    }
}

