#include "pressurewall.h"

#include "solidonsolidreaction.h"

using namespace kMC;


PressureWall::PressureWall(SolidOnSolidSolver &mutexSolver,
                       const double E0,
                       const double sigma0,
                       const double r0):
    SolidOnSolidEvent("MovingWall", "h0", true, true),
    m_mutexSolver(mutexSolver),
    m_r0(r0),
    m_s0(sigma0),
    m_E0(E0)
{
    BADAss(E0, <, 0, "E0 should be negative.");
}

PressureWall::~PressureWall()
{

}

void PressureWall::initialize()
{
    m_localPressure.set_size(solver()->length(), solver()->width());
    m_localPressure.zeros();

    BADAssBool(isActive() && hasStarted());

    m_thetaPrev = 0;

    for (uint x = 0; x < solver()->length(); ++x)
    {
        for (uint y = 0; y < solver()->width(); ++y)
        {
            m_thetaPrev += exp(solver()->height(x, y)/m_r0);
        }
    }

    m_thetaPrev /= solver()->area();

    m_height = m_r0*std::log(solver()->area()*m_s0*m_thetaPrev/(-m_E0));

    solver()->initializeReactions();

    BADAssClose(pressureEnergySum(), m_E0, 1E-5);

}

void PressureWall::reset()
{

    //Check that everything is updated correctly from previous runs.
#ifndef NDEBUG
    for (uint x = 0; x < solver()->length(); ++x)
    {
        for (uint y = 0; y < solver()->width(); ++y)
        {
            BADAssClose(localPressureEvaluate(x, y), localPressure(x, y), 1E-3);
        }
    }
#endif

    //Calculate new height of wall to conserve total force.
    _rescaleHeight();

    _updatePressureRates();

    BADAssClose(pressureEnergySum(), m_E0, 1E-5);

    if ((cycle() + 1)%10000 == 0)
    {
        recalculateAllPressures();
    }

    cout << "derp: update reactions which changed neighbors" << endl;

}

void PressureWall::_rescaleHeight()
{

    double theta = 0;
    for (uint x = 0; x < solver()->length(); ++x)
    {
        for (uint y = 0; y < solver()->width(); ++y)
        {
            theta += exp(solver()->height(x, y)/m_r0);
        }
    }

    theta /= solver()->area();

    m_heightChange = m_r0*std::log(theta/m_thetaPrev);

    m_thetaPrev = theta;

    m_height += m_heightChange;

    setValue(m_height - dependency("height")->value());

}

void PressureWall::_updatePressureRates()
{

    double rateChange;

    double expFac = expSmallArg(-m_heightChange/m_r0);

    DiffusionDeposition *reaction;

    for (uint x = 0; x < solver()->length(); ++x)
    {
        for (uint y = 0; y < solver()->width(); ++y)
        {

            reaction = &solver()->reaction(x, y);

            cout << "DERP: IF NOT AFFECTED" << endl;
            if (true)
            {

                rateChange = expSmallArg(-solver()->alpha()*m_localPressure(x, y)*(expFac - 1));

                //For every affected particle we update only those who include the pressure term.
                //Vector is set up in initialize based on virtual reaction function isPressureAffected().

                double prevDiffRate = reaction->diffusionRate();
                reaction->setDiffusionRate(prevDiffRate*rateChange);
                reaction->changeRate(reaction->diffusionRate() + reaction->depositionRate());

                m_localPressure(x, y) *= expFac;
            }
            else
            {
                m_localPressure(x, y) = localPressureEvaluate(x, y);
            }

            BADAssClose(localPressureEvaluate(x, y), m_localPressure(x, y), 1E-3, "incorrect pressure update", [&] ()
            {
                BADAssSimpleDump(cycle(), x, y, localPressure(x, y), expFac, m_heightChange);
            });
        }
    }
}

void PressureWall::recalculateAllPressures()
{
    for (uint x = 0; x < solver()->length(); ++x)
    {
        for (uint y = 0; y < solver()->width(); ++y)
        {
            m_localPressure(x, y) = localPressureEvaluate(x, y);
        }
    }
}

