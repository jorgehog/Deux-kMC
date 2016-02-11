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
    solver.registerObserver(this);
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
            const uint &x1 = csc.x1;
            const uint &y1 = csc.y1;

            m_potentialValues(x1, y1) = potentialFunction(x1, y1);
        }

        m_potentialValues(x, y) = potentialFunction(x, y);
    }

    else
    {
        const CurrentConfinementChange &ccc = solver().confiningSurfaceEvent().currentConfinementChange();

        double dh = solver().confiningSurfaceEvent().height() - ccc.prevHeight;

        m_expFac = expSmallArg(-dh/m_lD);

        for (uint x = 0; x < solver().length(); ++x)
        {
            for (uint y = 0; y < solver().width(); ++y)
            {
                SurfaceReaction &reaction = solver().surfaceReaction(x, y);

                double rateChange = expSmallArg(-solver().alpha()*m_potentialValues(x, y)*(m_expFac - 1));

                //For every affected particle we update only those who include the pressure term.
                //Vector is set up in initialize based on virtual reaction function isPressureAffected().

                double prevDiffRate = reaction.escapeRate();
                reaction.setEscapeRate(prevDiffRate*rateChange);

                m_potentialValues(x, y) *= m_expFac;

                BADAssClose(potentialFunction(x, y), m_potentialValues(x, y), 1E-5, "incorrect pressure update", [&] ()
                {
                    BADAssSimpleDump(x, y, rateChange, m_potentialValues(x, y), m_expFac, dh);
                });
            }
        }

    }
}

double RDLPotential::potentialFunction(const uint x, const uint y) const
{
    const double &h = solver().confiningSurfaceEvent().height();
    const int &hi = solver().height(x, y);

    return -m_s0*std::exp(-(h-hi-m_shift)/m_lD);
}
