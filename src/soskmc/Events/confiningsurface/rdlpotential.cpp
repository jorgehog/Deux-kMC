#include "rdlpotential.h"
#include "rdlsurface.h"
#include "sossolver.h"
#include "dissolutiondeposition.h"

RDLPotential::RDLPotential(SOSSolver &solver,
                           const double s0,
                           const double lD) :
    LocalCachedPotential(solver),
    m_s0(s0),
    m_lD(lD)
{

}

double RDLPotential::rdlEnergy(const double dh) const
{
    if (dh == 1)
    {
        return 0.;
    }

    return -m_s0*std::exp(-(dh-1)/m_lD);
}

double RDLPotential::rdlEnergyDeriv(const double dh) const
{
    return -rdlEnergy(dh)/m_lD;
}

double RDLPotential::expSmallArg(double arg)
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

void RDLPotential::notifyObserver(const kMC::Subjects &subject)
{
    if (subject == Subjects::SOLVER)
    {
        const CurrentSurfaceChange &csc = solver().currentSurfaceChange();

        const uint &x = csc.x;
        const uint &y = csc.y;

        if (csc.type == ChangeTypes::Double)
        {
            const int &x1 = csc.x1;
            const int &y1 = csc.y1;

            if (!solver().isOutsideBox(x1, y1))
            {
                m_potentialValues(x1, y1) = potentialFunction(x1, y1);
            }
        }

        m_potentialValues(x, y) = potentialFunction(x, y);
    }

    else
    {
        const CurrentConfinementChange &ccc = solver().confiningSurfaceEvent().currentConfinementChange();

        double dh = solver().confiningSurfaceEvent().height() - ccc.prevHeight;

        if (fabs(dh) < 1E-16)
        {
            return;
        }

        m_expFac = expSmallArg(-dh/m_lD);

        for (uint x = 0; x < solver().length(); ++x)
        {
            for (uint y = 0; y < solver().width(); ++y)
            {
                SurfaceReaction &reaction = solver().surfaceReaction(x, y);
                const double dhi = solver().confiningSurfaceEvent().height() - solver().height(x, y);

                if (m_potentialValues(x, y) == 0 || dhi == 1)
                {
                    m_potentialValues(x, y) = potentialFunction(x, y);
                    solver().registerAffectedReaction(&reaction);
                }

                else
                {
                    double rateChange = expSmallArg(-solver().alpha()*m_potentialValues(x, y)*(m_expFac - 1));

                    //For every affected particle we update only those who include the pressure term.
                    //Vector is set up in initialize based on virtual reaction function isPressureAffected().

                    double prevDiffRate = reaction.escapeRate();
                    reaction.setEscapeRate(prevDiffRate*rateChange);

                    m_potentialValues(x, y) *= m_expFac;

                    BADAssClose(potentialFunction(x, y), m_potentialValues(x, y), 1E-5, "incorrect pressure update", [&] ()
                    {
                        BADAssSimpleDump(x, y, dhi, rateChange, m_potentialValues(x, y), m_expFac, dh);
                    });
                }
            }
        }

    }
}

double RDLPotential::potentialFunction(const uint x, const uint y) const
{
    const double &h = solver().confiningSurfaceEvent().height();
    const int &hi = solver().height(x, y);

    return rdlEnergy(h - hi);
}
