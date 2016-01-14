#include "astarfirstpassage.h"

#include "../../sossolver.h"
#include "dissolutiondeposition.h"

#include "astarcpp/astarcpp.h"

#include "Events/confiningsurface/confiningsurface.h"

using namespace Tests;


AStarFirstPassage::AStarFirstPassage(SOSSolver &solver,
                                     const double maxdt,
                                     const int depositionBoxHalfSize,
                                     const double c) :
    Diffusion(solver, "AStarFirstPassage"),
    FirstPassageContinuum(solver, maxdt, depositionBoxHalfSize, c),
    m_boxSize(depositionBoxHalfSize*2 + 1),
    m_world(new World(m_boxSize, m_boxSize, m_boxSize)),
    m_pathFinder(new PathFinder(*m_world)),
    m_pathFindingJazzes(m_boxSize*m_boxSize),
    m_rootsOfInts(3)
{
    for (int i = 0; i < m_boxSize*m_boxSize; ++i)
    {
        m_pathFindingJazzes[i] = new PathFindingJazz;
    }

    m_rootsOfInts(0) = 1;
    m_rootsOfInts(1) = sqrt(2);
    m_rootsOfInts(2) = sqrt(3);
}


AStarFirstPassage::~AStarFirstPassage()
{
    delete m_world;
    delete m_pathFinder;

    for (PathFindingJazz *pfj : m_pathFindingJazzes)
    {
        delete pfj;
    }
}


void AStarFirstPassage::calculateLocalRatesAndUpdateDepositionRates()
{
    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            solver().surfaceReaction(x, y).setDepositionRate(0);
        }
    }

    const int &l = depositionBoxHalfSize();

    int xTrans;
    int yTrans;

    double r;

    SurfaceReaction *reaction;
    double localRate;

    const int hmax = solver().heights().max();

    const double D = DScaled();

    //
    bool dump = false;
    lammpswriter writer(4, "astar", "/tmp");
    lammpswriter surf(3, "surf", "/tmp");
    //

    if (dump)
    {
        surf.setSystemSize(2*l+1, 2*l+1, 2*l+1);
        writer.setSystemSize(2*l+1, 2*l+1, 2*l+1);
    }

    for (uint n = 0; n < nOfflatticeParticles(); ++n)
    {
        m_nPathFinds = 0;
        m_world->ResetBlocks();

        for (uint x = 0; x < solver().length(); ++x)
        {
            for (uint y = 0; y < solver().width(); ++y)
            {
                m_localRates(x, y, n) = 0;
            }
        }

        const double &zp = particlePositions(2, n);

        //if the particle is outside the pathfinding box of the
        //highest surface site it is not involved in any pathings.
        if (zp > hmax + l + 1)
        {
            continue;
        }

        const double &xp = particlePositions(0, n);
        const double &yp = particlePositions(1, n);

        const int ix = int(round(xp));
        const int iy = int(round(yp));
        const int iz = int(round(zp));

        const double dx = ix - xp;
        const double dy = iy - yp;
        const double dz = iz - zp;

        const double x0w = (l-dx);
        const double y0w = (l-dy);
        const double z0w = (l-dz);

        //Should find shortest path between xp yp zp and all heights inside a 2l+1 cube
        //Setup the cube first?
        //Yes.. then find all the paths and keep only those who have a total path length < l

        if (dump)
        {
            writer.initializeNewFile(n);
            surf.initializeNewFile(n);
            writer << 0 << l << l << l;
        }

        for (int xscan = -l; xscan <= l; ++xscan)
        {
            for (int yscan = -l; yscan <= l; ++yscan)
            {
                solver().boundaryLatticeTransform(xTrans, yTrans, ix + xscan, iy + yscan, iz);

                const int &h = solver().height(xTrans, yTrans);

                //if there is no room to deposit.
                if (h > solver().confiningSurfaceEvent().height() - 1)
                {
                    continue;
                }

                const double dz2 = (h+1-zp)*(h+1-zp);

                //if the surface site is outside the pathfinding box we skip it.
                if (dz2 < l*l)
                {
                    const int xw = xscan + l;
                    const int yw = yscan + l;
                    const int zw = h - iz + l;

                    for (int zwIter = 0; zwIter <= zw; ++zwIter)
                    {
                        m_world->MarkPosition(xw, yw, zwIter, true);

                        if (dump)
                        {
                            surf << xw << yw << zwIter;
                        }
                    }

//                    //we mark the site above as blocked as well since
//                    //we do not allow particles to deposit along the way
//                    //to the site. When pathfinding to site x,y,z will
//                    //be executed, we will unblock the target site.
//                    m_pathFinder->markBlockedPosition(xw, yw, zw+1);


                    m_pathFindingJazzes[m_nPathFinds]->xTrans = xTrans;
                    m_pathFindingJazzes[m_nPathFinds]->yTrans = yTrans;
                    m_pathFindingJazzes[m_nPathFinds]->xEnd = xw;
                    m_pathFindingJazzes[m_nPathFinds]->yEnd = yw;
                    m_pathFindingJazzes[m_nPathFinds]->zEnd = zw+1;
                    m_nPathFinds++;
                }
            }
        }

        for (uint i = 0; i < m_nPathFinds; ++i)
        {
            PathFindingJazz *pfj = m_pathFindingJazzes[i];

            xTrans = pfj->xTrans;
            yTrans = pfj->yTrans;

            int &xw = pfj->xEnd;
            int &yw = pfj->yEnd;
            int &zw = pfj->zEnd;

            r = 0;
            const SearchNode *crumb = m_pathFinder->findPath(l, l, l, xw, yw, zw);

            //no solution
            if (crumb == nullptr)
            {
                if (dump)
                {
                    writer << 3+i << l << l << l;
                }

                continue;
            }

            //We skip the start position
            //and calculate distance between the first and second node
            //explicitly
            crumb = crumb->next;

            const double r2FirstX = (crumb->x-x0w)*(crumb->x-x0w);
            const double r2FirstY = (crumb->y-y0w)*(crumb->y-y0w);
            const double r2FirstZ = (crumb->z-z0w)*(crumb->z-z0w);

            r += sqrt(r2FirstX + r2FirstY + r2FirstZ);

            while (crumb->next != nullptr)
            {
                double pathR2 = crumb->GetDistanceSquared(crumb->next->x,
                                                          crumb->next->y,
                                                          crumb->next->z);
                double pathR = m_rootsOfInts(pathR2 - 1);

                //so we don't calculate sqrt2 and sqrt3 over and over
                r += pathR;

                BADAssClose(pathR, sqrt(pathR2), 1E-10);

                if (dump)
                {
                    writer << 3+i << crumb->x << crumb->y << crumb->z;
                }

                crumb = crumb->next;
            }

            if (dump)
            {
                writer << 3+i << crumb->x << crumb->y << crumb->z;
            }


            reaction = &solver().surfaceReaction(xTrans, yTrans);
            localRate = localRateOverD(r*r);
            m_localRates(xTrans, yTrans, n) = localRate;
            reaction->setDepositionRate(reaction->depositionRate() + localRate*D);
        }

        if (dump)
        {

            if (surf.valueCounter() == 0)
            {
                surf << 0 << 0 << 0;
            }

            writer.finalize();
            surf.finalize();

        }
    }
}


void AStarFirstPassage::execute()
{
}
