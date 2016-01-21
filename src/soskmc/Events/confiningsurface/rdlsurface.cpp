#include "rdlsurface.h"

#include "../../dissolutiondeposition.h"
#include "../../concentrationboundaryreaction.h"

RDLSurface::RDLSurface(SOSSolver &solver,
                       const double E0,
                       const double s0,
                       const double lD) :
    ConfiningSurface(solver, "RDLSurface", "l0", true, true),
    m_lD(lD),
    m_s0(s0),
    m_E0(E0),
    m_RDLEnergy(solver.length(), solver.width(), fill::zeros),
    m_ratioPartialSums(solver.length(), solver.width())
{

}

RDLSurface::~RDLSurface()
{

}

double RDLSurface::expSmallArg(double arg)
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

double RDLSurface::partialThetaRatio(const uint x, const uint y) const
{
    return std::exp((solver().height(x, y) - m_r0LogThetaPrev)/m_lD);
}

double RDLSurface::RDLEnergySum() const
{
    double s = 0;
    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            BADAssClose(m_RDLEnergy(x, y), evaluateRDLEnergy(x, y), 1E-5);

            s += m_RDLEnergy(x, y);
        }
    }

    return s;
}

void RDLSurface::findNewHeight()
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

    setHeight(height() + m_heightChange);

    m_expFac = expSmallArg(-m_heightChange/m_lD);

    m_ratioPartialSums*=m_expFac;
}

void RDLSurface::updateRatesFor(SurfaceReaction &reaction)
{
    const uint x = reaction.x();
    const uint y = reaction.y();

    double rateChange = expSmallArg(-solver().alpha()*m_RDLEnergy(x, y)*(m_expFac - 1));

    //For every affected particle we update only those who include the pressure term.
    //Vector is set up in initialize based on virtual reaction function isPressureAffected().

    double prevDiffRate = reaction.escapeRate();
    reaction.setEscapeRate(prevDiffRate*rateChange);

    m_RDLEnergy(x, y) *= m_expFac;

    BADAssClose(evaluateRDLEnergy(x, y), m_RDLEnergy(x, y), 1E-5, "incorrect pressure update", [&] ()
    {
        BADAssSimpleDump(cycle(), x, y, rateChange, RDLEnergy(x, y), m_expFac, m_heightChange);
    });

//    BADAssClose(reaction.escapeRate(), reaction.calculateEscapeRate(), 1E-5, "incorrect rate update", [&] ()
//    {
//        double lp = RDLEnergy(x, y);
//        int h = solver().height(x, y);
//        BADAssSimpleDump(cycle(), x, y, rateChange, prevDiffRate, lp, m_expFac, m_heightChange, h, height());
//    });

}

double RDLSurface::bruteForceThetaRatio() const
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

void RDLSurface::setupTheta()
{
    double thetaShift = solver().averageHeight();
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

double RDLSurface::calculateKZrel(const double x, const double y, const double z, double &K, double &zRel) const
{
    const uint X = (uint)x;
    const uint Y = (uint)y;

    const int &h = solver().height(X, Y);
    const double &hl = height();

    const double D = hl - h;

    K = datum::pi/D;
    zRel = ((z-h) - D/2);

    return K*zRel;
}

void RDLSurface::_validateStoredEnergies() const
{
    //Check that everything is updated correctly from previous runs.
#ifndef NDEBUG
    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            BADAssClose(evaluateRDLEnergy(x, y), RDLEnergy(x, y), 1E-5);
        }
    }
#endif
}

void RDLSurface::recalculateAllRDLEnergies()
{
    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            recalculateRDLEnergy(x, y);
        }
    }
}

double RDLSurface::evaluateRDLEnergy(const uint x, const uint y) const
{
    return _RDLEnergyExpression(height() - solver().height(x, y));
}

void RDLSurface::execute()
{
    double value = height() - solver().averageHeight();

    setValue(value);

    if ((cycle() + 1)%100 == 0)
    {
        recalculateAllRDLEnergies();
    }
}

void RDLSurface::initialize()
{

}

void RDLSurface::reset()
{
    _validateStoredEnergies();

    BADAssClose(RDLEnergySum(), -m_E0, 1E-3);
}

void RDLSurface::initializeObserver(const Subjects &subject)
{
    (void) subject;

    setupTheta();

    setHeight(m_r0LogThetaPrev + m_lD*std::log(m_s0/m_E0));

    recalculateAllRDLEnergies();

    BADAssClose(RDLEnergySum(), -m_E0, 1E-5);
}

void RDLSurface::notifyObserver(const Subjects &subject)
{
    (void) subject;

    if (!hasStarted())
    {
        return;
    }

    const CurrentSurfaceChange &csc = solver().currentSurfaceChange();

    const uint &x = csc.x;
    const uint &y = csc.y;

    m_ratioPartialSums(x, y) = partialThetaRatio(x, y);
    m_RDLEnergy(x, y) = evaluateRDLEnergy(x, y);

    if (csc.type == ChangeTypes::Double)
    {
        const uint &x1 = csc.x1;
        const uint &y1 = csc.y1;

        m_ratioPartialSums(x1, y1) = partialThetaRatio(x1, y1);
        m_RDLEnergy(x1, y1) = evaluateRDLEnergy(x1, y1);
    }

    findNewHeight();

    for (uint _x = 0; _x < solver().length(); ++_x)
    {
        for (uint _y = 0; _y < solver().width(); ++_y)
        {
            if (_x == x && y == _y)
            {
                continue;
            }

            else if (csc.type == ChangeTypes::Double)
            {
                if (_x == csc.x1 && _y == csc.y1)
                {
                    continue;
                }
            }

            SurfaceReaction &reaction = mutexSolver().surfaceReaction(_x, _y);
            updateRatesFor(reaction);
        }
    }
}

bool RDLSurface::acceptDiffusionMove(const double x0, const double y0, const double z0,
                                     const double x1, const double y1, const double z1) const
{
    double kZrel0 = calculateKZrel(x0, y0, z0);
    double kZrel1 = calculateKZrel(x1, y1, z1);

    double a = pow(cos(kZrel0)/cos(kZrel1), 2.0);

    return rng.uniform() <= a;
}

double RDLSurface::diffusionDrift(const double x, const double y, const double z) const
{
    return 0;

    double K;
    double KZrel = calculateKZrel(x, y, z, K);

    return 2*K*tan(KZrel);
}

