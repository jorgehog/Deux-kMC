#include "rdlsurface.h"

#include "../../dissolutiondeposition.h"
#include "../../concentrationboundaryreaction.h"
#include "rdlpotential.h"

RDLSurface::RDLSurface(SOSSolver &solver,
                       RDLPotential &potential,
                       const double Pl) :
    ConfiningSurface(solver, "RDLSurface", "l0", true, true),
    m_potential(potential),
    m_Pl(Pl),
    m_ratioPartialSums(solver.length(), solver.width())
{
    registerObserver(&potential);
}

RDLSurface::~RDLSurface()
{

}

double RDLSurface::partialThetaRatio(const uint x, const uint y) const
{
    return std::exp((solver().height(x, y) - m_ldLogThetaPrev)/m_potential.lD());
}

double RDLSurface::RDLEnergySum() const
{
    double s = 0;
    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            s += m_potential.potentialValue(x, y);
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
                BADAssSimpleDump(m_ratioPartialSums(x, y), partialThetaRatio(x, y), m_ldLogThetaPrev);
            });
        }
    }

    const double heightChange = m_potential.lD()*std::log(ratio);

    BADAss(abs(heightChange), <=, 1, "undefined behavior when wall is moved more than one cell per iteration.");

    m_ldLogThetaPrev += heightChange;

    setHeight(height() + heightChange);

    //this replaces old m_r0LogThetaPrev with the new value.
    m_ratioPartialSums *= m_potential.expFac();
}

double RDLSurface::bruteForceThetaRatio() const
{
    double ratio = 0;

    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            ratio += std::exp(((long double)solver().height(x, y) - m_ldLogThetaPrev)/m_potential.lD());
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
            thetaPrev += std::exp((solver().height(x, y) - thetaShift)/m_potential.lD());
        }
    }

    m_ldLogThetaPrev = m_potential.lD()*log(thetaPrev) + thetaShift;


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

void RDLSurface::execute()
{
    setValue(height() - solver().averageHeight());
}

void RDLSurface::initialize()
{

}

void RDLSurface::reset()
{
    BADAssClose(RDLEnergySum(), -m_Pl*solver().area(), 1E-3);
}

void RDLSurface::initializeObserver(const Subjects &subject)
{
    (void) subject;

    setupTheta();

    setHeight(RDLPotential::m_shift + m_ldLogThetaPrev + m_potential.lD()*std::log(m_potential.s0()/(m_Pl*solver().area())));

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

    if (csc.type == ChangeTypes::Double)
    {
        const uint &x1 = csc.x1;
        const uint &y1 = csc.y1;

        m_ratioPartialSums(x1, y1) = partialThetaRatio(x1, y1);
    }

    findNewHeight();

    BADAssClose(RDLEnergySum(), -m_Pl*solver().area(), 1E-5);
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

