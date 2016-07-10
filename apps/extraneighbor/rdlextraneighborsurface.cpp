#include "rdlextraneighborsurface.h"


RDLExtraNeighborSurface::RDLExtraNeighborSurface(SOSSolver &solver,
                                                 RDLPotential &rdlPotential,
                                                 ExtraNeighbor &extraNeighborPotential,
                                                 const double Pl) :
    ConfiningSurface(solver, "RDLExtraNeighborSurface"),
    m_rdlPotential(rdlPotential),
    m_extraNeighborPotential(extraNeighborPotential),
    m_Pl(Pl),
    m_nFactor(extraNeighborPotential.relBondEnergy()*rdlPotential.lD()*(1 - exp(-1/rdlPotential.lD()))/20),
    m_prevMinSet(false),
    m_ldExpFac(exp(-1./rdlPotential.lD()))
{
    registerObserver(&rdlPotential);
    registerObserver(&extraNeighborPotential);
}

double RDLExtraNeighborSurface::getRdlEquilibrium() const
{
    const double &ld = m_rdlPotential.lD();
    const double &s0 = m_rdlPotential.s0();

    long double m;
    if (solver().initialized())
    {
        m = solver().averageHeight();
    }

    else
    {
        m = arma::accu(solver().heights())/solver().area();
    }

    long double T = 0;
    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            T += exp((solver().height(x, y) - m)/ld);
        }
    }

    T /= solver().area();

    return 1 + m + ld*log(s0*T/m_Pl);

}

void RDLExtraNeighborSurface::dumpProfile(const uint n) const
{
    const int m = solver().heights().max();

    vec hls = linspace(m+1, m+4);

    mat res(2, hls.size());

    for (uint i = 0; i < hls.size(); ++i)
    {
        res(0, i) = hls(i);
        res(1, i) = totalForce(hls(i));
    }

    string ending = to_string(n) + ".arma";

    res.save("/tmp/mechEq" + ending);
}

void RDLExtraNeighborSurface::findNewHeight()
{
    const double h = getHeightBisection();

    BADAss(h, >=, solver().heights().max() + 1 - 1E-10, "hmm", [&] {
        const double hb = getHeightBisection();
        cout << hb << endl;
    });

    setHeight(h);
}

double RDLExtraNeighborSurface::totalForce(const double hl) const
{
    double rdlContribution = 0;
    double extraNeighborContribution = 0;

    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            const int &hi = solver().height(x, y);

            const double dh = hl - hi;
            rdlContribution += m_rdlPotential.rdlEnergy(dh);

            if (dh < 2)
            {
                const double dh2 = dh*dh;
                extraNeighborContribution += 128/(dh2*dh2*dh2*dh) - 1;
            }
        }
    }

    extraNeighborContribution /= solver().area(); //we devide by area since we have devided the entire force balance equation by the area since we use sigma0 and not z0
    rdlContribution /= solver().area();

    return m_Pl + rdlContribution + extraNeighborContribution*m_nFactor;
}

double RDLExtraNeighborSurface::totalForceDeriv(const double hl) const
{
    double rdlContribution = 0;
    double extraNeighborContribution = 0;

    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            const int &hi = solver().height(x, y);

            const double dh = hl - hi;
            rdlContribution += m_rdlPotential.rdlEnergyDeriv(dh);

            if (dh < 2)
            {
                const double dh2 = dh*dh;
                extraNeighborContribution += -896/(dh2*dh2*dh2*dh2);
            }
        }
    }

    extraNeighborContribution /= solver().area();
    rdlContribution /= solver().area();

    return rdlContribution + extraNeighborContribution*m_nFactor;
}

//! First we find a negative point,
//! and then we find the far equilibrium by bisection
//! then we decide which equilibrium to choose
//! based on the current force.
double RDLExtraNeighborSurface::getHeightBisection()
{
    bool shouldRepel = false;

    const double &currentHeight = height();

    //this is the minimum height the solution can have
    const int contactHeight = solver().heights().max() + 1;

    const double rdlEquilibrium = getRdlEquilibrium();

    if (fabs(contactHeight - currentHeight) < 1E-10)
    {
        if (rdlEquilibrium < contactHeight + 1)
        {
            return contactHeight;
        }

        double repcurrent = 0;
        double rep_energy_barrier = 0;
        double e_bonds = 0;

        for (uint x = 0; x < solver().length(); ++x)
        {
            for (uint y = 0; y < solver().width(); ++y)
            {
                const int &hi = solver().height(x, y);
                if (contactHeight - hi == 1)
                {
                    rep_energy_barrier += m_rdlPotential.s0();
                    e_bonds += 1;
                }

                else
                {
                    rep_energy_barrier -= m_rdlPotential.potential(x, y);
                }

                repcurrent -= m_rdlPotential.potential(x, y);
            }
        }

        rep_energy_barrier *= m_ldExpFac;
        const double xi = 1 - m_ldExpFac;

        const double wF = m_Pl/(m_rdlPotential.lD())*solver().area();

        if (rep_energy_barrier - repcurrent > wF + e_bonds*xi)
        {
            shouldRepel = true;
        }

        else
        {
            return currentHeight;
        }
    }

    else
    {

        //if we have a net repulsion, we always choose the far point.
        const double currentForce = totalForce(currentHeight);

        //we are already in equilibrium
        if (fabs(currentForce) < 1E-10)
        {
            return currentHeight;
        }

        shouldRepel = currentForce < 0;

        if ((rdlEquilibrium > contactHeight + 1) && shouldRepel)
        {
            return rdlEquilibrium;
        }

    }

    //-- Finding a point with negative force (repulsion)

    const double contactForce = totalForce(contactHeight+0.01);

    //if the maximum left value is negative there
    //is no need to use bisection to obtain a negative value
    double hNegative;
    double minForce;
    if (contactForce < 0)
    {
        hNegative = contactHeight;
        minForce = contactForce;
    }

    else
    {
        bool bisectForMin = true;
        if (m_prevMinSet)
        {
            const double prevForce = totalForce(m_prevMin);
            if (prevForce < 0)
            {
                hNegative = m_prevMin;
                bisectForMin = false;
                minForce = prevForce;
            }
        }

        if (bisectForMin)
        {
            hNegative = bisectMinima(contactHeight,
                                     contactHeight + 2,
                                     totalForceDeriv(contactHeight),
                                     1E-3);

            m_prevMinSet = true;
            m_prevMin = hNegative;
            minForce = totalForce(hNegative);

            //no equilibriums (roots)
            if (minForce > 0)
            {
                return contactHeight;
            }
        }
    }

    //--

    //--Finding the far height

    double farHeight;

    //outside range of interactions
    if (hNegative > contactHeight + 1)
    {
        farHeight = rdlEquilibrium;
    }

    else
    {
        //DERP
        double hPositive = hNegative+1;
        while (totalForce(hPositive) < 0)
        {
            hPositive++;
        }

        farHeight = bisect(hNegative, hPositive, minForce);
    }

    //--

    //repelling to the far height
    if (shouldRepel)
    {
        return farHeight;
    }

    //if we get here we are guaranteed an attraction
    else
    {
        //in this case we are attracted untill we reach the far point
        if ((currentHeight > farHeight) && (farHeight > contactHeight))
        {
            return farHeight;
        }

        //only option left is to be attracted to the surface
        else
        {
            return contactHeight;
        }
    }
}

double RDLExtraNeighborSurface::bisect(double min,
                                       double max,
                                       double fmin,
                                       const double eps,
                                       const uint nMax) const
{
    BADAss(min, <, max, "hmm", [&] ()
    {
        dumpProfile(0);
    });

    double fmid;
    double mid = 0;

    uint n = 0;
    while ((n < nMax) && (max - min > eps))
    {
        mid = 0.5*(min + max);
        fmid = totalForce(mid);

        if ((mid == max) || (mid == min) || fmid == 0)
        {
            break;
        }

        else if (fmid*fmin < 0)
        {
            max = mid;
        }

        else
        {
            min = mid;
            fmin = fmid;
        }

        n++;
    }

    return mid;
}


double RDLExtraNeighborSurface::bisectMinima(double min,
                                             double max,
                                             double fmin,
                                             const double eps,
                                             const uint nMax) const
{
    BADAss(min, <, max, "hmm", [&] ()
    {
        dumpProfile(0);
    });

    double fmid;
    double mid = 0;

    uint n = 0;
    while ((n < nMax) && (max - min > eps))
    {
        mid = 0.5*(min + max);
        fmid = totalForceDeriv(mid);

        if ((mid == max) || (mid == min) || fmid == 0)
        {
            break;
        }

        else if (fmid*fmin < 0)
        {
            max = mid;
        }

        else
        {
            min = mid;
            fmin = fmid;
        }

        n++;
    }

    return mid;
}

void RDLExtraNeighborSurface::initializeObserver(const Subjects &subject)
{
    (void) subject;

    m_prevMinSet = false;

    findNewHeight();
}

void RDLExtraNeighborSurface::notifyObserver(const Subjects &subject)
{
    (void) subject;

    if (solver().currentSurfaceChange().type == ChangeTypes::Single)
    {
        findNewHeight();
    }

    //we perform checks if the surfaces are not in contact if the forces are indeed balanced.
    if (fabs(height() - (solver().heights().max() + 1) > 1E-10))
    {
        BADAssClose(0, totalForce(height()), 1E-10, "hmm", [&] ()
        {
            BADAssSimpleDump(solver().heights().max(), height());
            dumpProfile(0);
        });
    }
}

void RDLExtraNeighborSurface::execute()
{

}

bool RDLExtraNeighborSurface::hasSurface() const
{
    return true;
}

bool RDLExtraNeighborSurface::acceptDiffusionMove(const double x0, const double y0, const double z0, const double x1, const double y1, const double z1) const
{
    (void) x0;
    (void) y0;
    (void) z0;
    (void) x1;
    (void) y1;
    (void) z1;

    return true;
}

double RDLExtraNeighborSurface::diffusionDrift(const double x, const double y, const double z) const
{
    (void) x;
    (void) y;
    (void) z;

    return 0.0;
}
