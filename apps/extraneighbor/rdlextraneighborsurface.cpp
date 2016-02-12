#include "rdlextraneighborsurface.h"

RDLExtraNeighborSurface::RDLExtraNeighborSurface(SOSSolver &solver,
                                                 RDLPotential &rdlPotential,
                                                 ExtraNeighbor &extraNeighborPotential,
                                                 const double Pl) :
    ConfiningSurface(solver, "RDLExtraNeighborSurface"),
    m_rdlPotential(rdlPotential),
    m_extraNeighborPotential(extraNeighborPotential),
    m_Pl(Pl)
{
    registerObserver(&rdlPotential);
    registerObserver(&extraNeighborPotential);
}

void RDLExtraNeighborSurface::findNewHeight()
{
    const double eps = 1E-10;

    double h = height();
    const int hMax = solver().heights().max();

    //standard RDL height since the n-potential is zero.
    if (h > hMax + 2)
    {
        setHeight(iteratingExpression(h));
        return;
    }

    double hPrev = h + 2*eps;

//    cout << h << endl;

    vector<double> savedHeights;
    uint i = 1;

    while (fabs(h - hPrev) > eps)
    {
        hPrev = h;

        h = iteratingExpression(hPrev);

//        cout << h << endl;

        //checks the new time step against previous time steps
        //if it is identical, then the iteration failed.
        for (const double &savedHeight : savedHeights)
        {
            if (fabs(savedHeight - h) < 1E-15)
            {
//                cout << "Fixed point iteration failed." << savedHeight << " " << h << endl;

//                for (const double &savedHeight2 : savedHeights)
//                {
//                    cout << savedHeight2 << " ";
//                }
//                cout << endl;
//                sleep(3);
                setHeight(hMax + 1);
                return;
            }
        }


        //cycle previous time steps back one index and add the new one to the end.
        if (!savedHeights.empty())
        {
            for (uint j = 0; j < savedHeights.size() - 1; ++j)
            {
                savedHeights.at(j) = savedHeights.at(j + 1);
            }
            savedHeights.back() = h;
        }

        //every 100 iteration cycles we check against one extra time step
        //to be able to locate larger loops, i.e. a-b-c-d-e-f-d-e-f-d-e-f-.. etc.
        //and not only let's say a-b-c-b-c-b-c..
        if (i % 20 == 0)
        {
            savedHeights.push_back(h);
        }

        i++;

    }

//    cout << "fin" << endl;

    //if the eq distance is below the max surface height,
    //we set the surface to rest on top of one another.
    h = h < (hMax + 1) ? (hMax + 1) : h;

    setHeight(h);
}

double RDLExtraNeighborSurface::iteratingExpression(const double hl) const
{
    const double &lD = m_rdlPotential.lD();
    const double &s0 = m_rdlPotential.s0();

    const double nFactor = lD*(s0 + (1-exp(-1/lD))/solver().area())/20;

    double nSum = 0;
    double theta = 0;

    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            const int &hi = solver().height(x, y);

            theta += exp(hi/m_rdlPotential.lD());

            if (hl - hi < 2 && hl >= hi + 1)
            {
                nSum += 128/pow(hl - hi, 7.0) - 1;
            }

        }
    }

    theta /= solver().area();

    return 1 - lD*log((m_Pl+nSum*nFactor)/(s0*theta));
}

double RDLExtraNeighborSurface::mechanicalEquilibriumCondition(const double hl) const
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
        }

    }


    return m_Pl - rdlContribution + extraNeighborContribution;
}

void RDLExtraNeighborSurface::initializeObserver(const Subjects &subject)
{
    (void) subject;

    setHeight(solver().heights().max() + 1);
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
