#include "rdlextraneighborsurface.h"


RDLExtraNeighborSurface::RDLExtraNeighborSurface(SOSSolver &solver,
                                                 RDLPotential &rdlPotential,
                                                 ExtraNeighbor &extraNeighborPotential,
                                                 const double Pl) :
    ConfiningSurface(solver, "RDLExtraNeighborSurface"),
    m_rdlPotential(rdlPotential),
    m_extraNeighborPotential(extraNeighborPotential),
    m_Pl(Pl),
    m_nFactor(rdlPotential.lD()*(1 - exp(-1/rdlPotential.lD()) + rdlPotential.s0())/20)
{
    registerObserver(&rdlPotential);
    registerObserver(&extraNeighborPotential);
}

double RDLExtraNeighborSurface::getRdlEquilibrium()
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

void RDLExtraNeighborSurface::dumpProfile() const
{
    vec hls = linspace(1, 5);
    vec Fs(hls.size());

    const int m = solver().heights().max();

    for (uint i = 0; i < hls.size(); ++i)
    {
        Fs(i) = totalForce(hls(i) + m);
    }

    hls.save("/tmp/mechEqHls.arma");
    Fs.save("/tmp/mechEqFs.arma");
}

void RDLExtraNeighborSurface::findNewHeight()
{
    bisection();

    //we perform checks if the surfaces are not in contact if the forces are indeed balanced.
    if (height() != solver().heights().max() + 1)
    {
        BADAssClose(0, totalForce(height()), 1E-10, "hmm", [&] ()
        {
            dumpProfile();
            BADAssSimpleDump(totalAttraction(), m_Pl);
        });
    }
}

//double RDLExtraNeighborSurface::iteratingExpression(const double hl) const
//{
//    const double &lD = m_rdlPotential.lD();
//    const double &s0 = m_rdlPotential.s0();

//    const double nFactor = lD*(s0 + (1-exp(-1/lD))/solver().area())/20;

//    double nSum = 0;
//    double theta = 0;

//    for (uint x = 0; x < solver().length(); ++x)
//    {
//        for (uint y = 0; y < solver().width(); ++y)
//        {
//            const int &hi = solver().height(x, y);

//            theta += exp(hi/m_rdlPotential.lD());

//            if (hl - hi < 2 && hl >= hi + 1)
//            {
//                nSum += 128/pow(hl - hi, 7.0) - 1;
//            }

//        }
//    }

//    theta /= solver().area();

//    return 1 - lD*log((m_Pl+nSum*nFactor)/(s0*theta));
//}

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

    extraNeighborContribution /= solver().area();
    rdlContribution /= solver().area();

    return m_Pl + rdlContribution + extraNeighborContribution*m_nFactor;
}

double RDLExtraNeighborSurface::totalRepulsion() const
{
    return m_rdlPotential.sum();
}

double RDLExtraNeighborSurface::totalAttraction() const
{
    double extraNeighborContribution = 0;

    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            const int &hi = solver().height(x, y);

            const double dh = height() - hi;

            if (dh < 2)
            {
                extraNeighborContribution += 128/pow(dh, 7) - 1;
            }
        }
    }

    extraNeighborContribution /= solver().area();

    return m_Pl + extraNeighborContribution*m_nFactor;
}

//void RDLExtraNeighborSurface::fixPointIteration()
//{
//    BADAssBreak("");

//    const double eps = 1E-10;

//    double h = height();
//    const int hMax = solver().heights().max();

//    //standard RDL height since the n-potential is zero.
//    if (h > hMax + 2)
//    {
//        setHeight(iteratingExpression(h));
//        return;
//    }

//    double hPrev = h + 2*eps;

//    vector<double> savedHeights;
//    uint i = 1;

//    while (fabs(h - hPrev) > eps)
//    {
//        hPrev = h;

//        h = iteratingExpression(hPrev);

//        //checks the new time step against previous time steps
//        //if it is identical, then the iteration failed.
//        for (const double &savedHeight : savedHeights)
//        {
//            if (fabs(savedHeight - h) < 1E-15)
//            {
//                setHeight(hMax + 1);
//                return;
//            }
//        }


//        //cycle previous time steps back one index and add the new one to the end.
//        if (!savedHeights.empty())
//        {
//            for (uint j = 0; j < savedHeights.size() - 1; ++j)
//            {
//                savedHeights.at(j) = savedHeights.at(j + 1);
//            }
//            savedHeights.back() = h;
//        }

//        //every 100 iteration cycles we check against one extra time step
//        //to be able to locate larger loops, i.e. a-b-c-d-e-f-d-e-f-d-e-f-.. etc.
//        //and not only let's say a-b-c-b-c-b-c..
//        if (i % 20 == 0)
//        {
//            savedHeights.push_back(h);
//        }

//        i++;

//    }

//    //if the eq distance is below the max surface height,
//    //we set the surface to rest on top of one another.
//    h = h < (hMax + 1) ? (hMax + 1) : h;

//    setHeight(h);
//}

void RDLExtraNeighborSurface::bisection()
{
    //this is the minimum height the solution can have
    const double contactHeight = solver().heights().max() + 1;

    //this is the height at which there are no particle-surface attraction
    const double noAttractionHeight = contactHeight + 1;

    //this is the minimum value of the force
    const double fNoAttraction = totalForce(noAttractionHeight);

    //In this case there are no equilibriums available above the contact
    if (fNoAttraction > 0)
    {
        setHeight(contactHeight);
        return;
    }

    double farHeight = getRdlEquilibrium();
    BADAssClose(bisect(noAttractionHeight,
                       noAttractionHeight + 1 + m_rdlPotential.lD()*log(m_rdlPotential.s0()/m_Pl),
                       fNoAttraction),
                farHeight,
                1E-5);

    const double &currentHeight = height();

    //if we have a net repulsion, we always choose the far point.
    if (totalForce(currentHeight) < 0)
    {
        setHeight(farHeight);
    }

    else
    {
        //in this case we are attracted untill we reach the far point
        if (currentHeight > farHeight)
        {
            setHeight(farHeight);
        }

        else
        {
            setHeight(contactHeight);
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
        dumpProfile();
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

void RDLExtraNeighborSurface::initializeObserver(const Subjects &subject)
{
    (void) subject;

    findNewHeight();
}

void RDLExtraNeighborSurface::notifyObserver(const Subjects &subject)
{
    (void) subject;

    findNewHeight();
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
