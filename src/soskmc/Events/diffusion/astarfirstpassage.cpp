#include "astarfirstpassage.h"

#include "../../sossolver.h"
#include "dissolutiondeposition.h"

#include "pathfinder.h"

AStarFirstPassage::AStarFirstPassage(SOSSolver &solver,
                                     const double maxdt,
                                     const int depositionBoxHalfSize,
                                     const double c) :
    Diffusion(solver, "AStarFirstPassage"),
    FirstPassageContinuum(solver, maxdt, depositionBoxHalfSize, c),
    m_boxSize(depositionBoxHalfSize*2 + 1),
    m_pathFinder(new PathFinder(m_boxSize, m_boxSize, m_boxSize)),
    m_pathFindingJazzes(m_boxSize*m_boxSize)
{
    for (int i = 0; i < m_boxSize*m_boxSize; ++i)
    {
        m_pathFindingJazzes[i] = new PathFindingJazz;
    }
}


AStarFirstPassage::~AStarFirstPassage()
{
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
            solver().surfaceReaction(x, y).setEscapeRate(0);
        }
    }

    const int &l = depositionBoxHalfSize();

    int xTrans;
    int yTrans;

    SurfaceReaction *r;
    double localRate;

    const int hmax = solver().heights().max();

    const double D = DScaled();

    //    lammpswriter writer(4, "astar", "/tmp");
    //    lammpswriter surf(3, "surf", "/tmp");
    //    surf.setSystemSize(2*l+1, 2*l+1, 2*l+1);
    //    writer.setSystemSize(2*l+1, 2*l+1, 2*l+1);

    std::vector<void *> path;
    float cost;

    for (uint n = 0; n < nOfflatticeParticles(); ++n)
    {
        m_nPathFinds = 0;

        for (uint x = 0; x < solver().length(); ++x)
        {
            for (uint y = 0; y < solver().width(); ++y)
            {
                m_localRates(x, y, n) = 0;
            }
        }

        const double &zp = particlePositions(2, n);

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
        //        writer.initializeNewFile(n);
        //        surf.initializeNewFile(n);
        //        writer << 0 << l << l << l;

        for (int xscan = -l; xscan <= l; ++xscan)
        {
            for (int yscan = -l; yscan <= l; ++yscan)
            {
                solver().boundaryLatticeTransform(xTrans, yTrans, ix + xscan, iy + yscan, iz);

                const int &h = solver().height(xTrans, yTrans);

                const double dz2 = (h+1-zp)*(h+1-zp);


                if (dz2 < l*l)
                {
                    const int xw = xscan + l;
                    const int yw = yscan + l;
                    const int zw = h - iz + l;

                    for (int zwIter = 0; zwIter <= zw; ++zwIter)
                    {
                        m_pathFinder->markBlockedPosition(xw, yw, zwIter);
                        //                        surf << xw << yw << zwIter;
                    }

                    m_pathFindingJazzes[m_nPathFinds]->xTrans = xTrans;
                    m_pathFindingJazzes[m_nPathFinds]->yTrans = yTrans;
                    m_pathFindingJazzes[m_nPathFinds]->xEnd = xw;
                    m_pathFindingJazzes[m_nPathFinds]->yEnd = yw;
                    m_pathFindingJazzes[m_nPathFinds]->zEnd = zw;
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

            int exitValue = m_pathFinder->solve(path, cost,
                                                l, l, l,
                                                xw, yw, zw);

            //            writer << 3+i << crumb->position->X << crumb->position->Y << crumb->position->Z;
            if (exitValue == micropather::MicroPather::NO_SOLUTION)
            {
                continue;
            }

            //We skip the start position
            //and calculate distance between the first and second node
            //explicitly using the double position starting point
            //x0w y0w z0w instead

            int x, y, z;
            m_pathFinder->getXYZ(x, y, z, path.at(1));

            const double r2FirstX = (x-x0w)*(x-x0w);
            const double r2FirstY = (y-y0w)*(y-y0w);
            const double r2FirstZ = (z-z0w)*(z-z0w);

            const double r2FirstXp = (x-l)*(x-l);
            const double r2FirstYp = (y-l)*(y-l);
            const double r2FirstZp = (z-l)*(z-l);

            const double rFirst = sqrt(r2FirstX + r2FirstY + r2FirstZ);
            const double rFirstP = sqrt(r2FirstXp + r2FirstYp + r2FirstZp);

            const double rSingle = cost - rFirstP + rFirst;
            const double r2 = rSingle*rSingle;

            r = &solver().surfaceReaction(xTrans, yTrans);
            localRate = localRateOverD(r2);
            m_localRates(xTrans, yTrans, n) = localRate;
            r->setEscapeRate(r->depositionRate() + localRate*D);
        }



        //        if (surf.valueCounter() == 0)
        //        {
        //            surf << 0 << 0 << 0;
        //        }

        //        writer.finalize();
        //        surf.finalize();
    }
}


void AStarFirstPassage::execute()
{
}
