#include "offlatticemontecarlo.h"

#include "../confiningsurface/confiningsurface.h"
#include "concentrationboundaryreaction.h"

#include "dissolutiondeposition.h"

#include "/home/jorgehog/code/astarcpp/astarcpp.h"

using namespace Tests;


OfflatticeMonteCarlo::OfflatticeMonteCarlo(SOSSolver &solver,
                                           const double maxdt,
                                           const int depositionBoxHalfSize,
                                           string type,
                                           string unit,
                                           bool hasOutput,
                                           bool storeValue) :
    Diffusion(solver, type, unit, hasOutput, storeValue),
    m_maxdt(maxdt),
    m_depositionBoxHalfSize(depositionBoxHalfSize),
    m_boxSize(depositionBoxHalfSize*2 + 1),
    m_world(new World(m_boxSize, m_boxSize, m_boxSize)),
    m_pathFinder(new PathFinder(*m_world)),
    m_pathFindingJazzes(m_boxSize*m_boxSize),
    m_particlePositions(3, 100),
    m_F(3, 100, fill::zeros),
    m_localRates(solver.length(), solver.width(), 100),
    m_localRatesForSite(100),
    m_nParticles(0)
{
    for (int i = 0; i < m_boxSize*m_boxSize; ++i)
    {
        m_pathFindingJazzes[i] = new PathFindingJazz;
    }
}

OfflatticeMonteCarlo::~OfflatticeMonteCarlo()
{
    delete m_world;
    delete m_pathFinder;

    for (PathFindingJazz *pfj : m_pathFindingJazzes)
    {
        delete pfj;
    }
}


void OfflatticeMonteCarlo::dump(const uint frameNumber, const string path) const
{
    Diffusion::dump(frameNumber, path);
    dumpDiffusingParticles(frameNumber, path);
}

uint OfflatticeMonteCarlo::dissolutionPaths(const uint x, const uint y) const
{
    bool lockedIn = solver().nNeighbors(x, y) == 5;
    bool closedTop = (solver().confiningSurfaceEvent().height() - solver().height(x, y)) <= 1;
    if ( lockedIn && closedTop )
    {
        return 0u;
    }

    return 1u;
}

void OfflatticeMonteCarlo::notifyObserver(const Subjects &subject)
{
    (void) subject;

    if (!hasStarted())
    {
        return;
    }

    const uint &x = solver().currentSurfaceChange().x;
    const uint &y = solver().currentSurfaceChange().y;
    const int &value = solver().currentSurfaceChange().value;

    //Remove particle based on the probability of it being the deposited
    if (value == 1)
    {
        double Rtot = 0;
        for (uint n = 0; n < nOfflatticeParticles(); ++n)
        {
            //use old rates to calculate probability of depositing
            m_localRatesForSite(n) = localRates(x, y, n);
            Rtot += m_localRatesForSite(n);
        }

        const uint N = chooseFromTotalRate(m_localRatesForSite.memptr(), nOfflatticeParticles(), Rtot);

        removeParticle(N);
    }

    //Add a particle on a unit sphere around the site;
    else
    {
        int dx, dy, dz;
        //        double x1, y1, z1;
        int x1, y1, z1;

        const int &z = solver().height(x, y) + 1;

        //        do
        //        {
        //            double theta = rng.uniform()*datum::pi;
        //            double phi = rng.uniform()*datum::pi*2;

        //            dx = sin(theta)*cos(phi);
        //            dy = sin(theta)*sin(phi);
        //            dz = cos(theta);

        //            x1 = solver().boundaryTransform(x, y, z, dx, 0);
        //            y1 = solver().boundaryTransform(x, y, z, dy, 1);
        //            z1 = z + dz;

        //            BADAssClose(1, dx*dx + dy*dy + dz*dz, 1E-3);


        //        } while (solver().isBlockedPosition(x1, y1, z1));

        solver().getRandomSolutionSite(x, y, z, dx, dy, dz);

        //        x1 = solver().boundaryTransform(x, y, z, dx, 0);
        //        y1 = solver().boundaryTransform(x, y, z, dy, 1);

        solver().boundaryLatticeTransform(x1, y1, x + dx, y + dy, z);

        z1 = z + dz;

        insertParticle(x1, y1, z1);
    }

    for (uint n = 0; n < nOfflatticeParticles(); ++n)
    {
        const double &x0 = particlePositions(0, n);
        const double &y0 = particlePositions(1, n);
        const double &z0 = particlePositions(2, n);

        //if the height change made the particle blocked,
        //we shift it a little and recalculate the rates
        if (solver().isBlockedPosition(x0, y0, z0))
        {
            uint dim;
            double delta;

            scanForDisplacement(n, dim, delta);
            m_particlePositions(dim, n) += delta;

            //            for (uint _x = 0; _x < solver().length(); ++_x)
            //            {
            //                for (uint _y = 0; _y < solver().width(); ++_y)
            //                {
            //                    m_localRates(_x, _y, n) = calculateLocalRateOverD(_x, _y, n);
            //                }
            //            }
        }

        else
        {
            //else we just recalulate the rate for the new height
            //            m_localRates(x, y, n) = calculateLocalRateOverD(x, y, n);
        }
    }

    //    DissolutionDeposition *r;
    //    for (uint _x = 0; _x < solver().length(); ++_x)
    //    {
    //        for (uint _y = 0; _y < solver().width(); ++_y)
    //        {
    //            r = &solver().surfaceReaction(_x, _y);

    //            const double newDepositionRate = depositionRate(_x, _y);

    //            r->setDepositionRate(newDepositionRate);
    //            r->changeRate(newDepositionRate + r->dissolutionRate());
    //        }
    //    }

}

double OfflatticeMonteCarlo::depositionRate(const uint x, const uint y) const
{
    double ROverD = 0;

    for (uint n = 0; n < nOfflatticeParticles(); ++n)
    {
        const double Rn = localRates(x, y, n);

        //        BADAssClose(Rn, calculateLocalRateOverD(x, y, n), 1E-3);

        ROverD += Rn;
    }

    //Dscaled ensures that this rate has the proper time unit
    return ROverD*DScaled();
}

void OfflatticeMonteCarlo::executeDiffusionReaction(SOSDiffusionReaction *reaction, const int x, const int y, const int z)
{
    (void) reaction;
    (void) x;
    (void) y;
    (void) z;

    throw std::logic_error("invalid diffusion model.");
}

void OfflatticeMonteCarlo::executeConcentrationBoundaryReaction(ConcentrationBoundaryReaction *reaction)
{
    const uint n =reaction->freeBoundarySites();
    const uint nChosen = rng.uniform()*n;

    uint ixi;
    int iz;
    reaction->getFreeBoundarSite(nChosen, ixi, iz);

    //place the particle randomly inside the selected box.
    double z = iz + (rng.uniform() - 0.5);
    double xi = ixi + (rng.uniform() - 0.5);

    if (reaction->dim() == 0)
    {
        insertParticle(reaction->location(), xi, z);
    }

    else
    {
        insertParticle(xi, reaction->location(), z);
    }
}

bool OfflatticeMonteCarlo::isBlockedPosition(const uint x, const uint y, const int z) const
{
    (void) x;
    (void) y;
    (void) z;

    //derp: hard core potential
    return false;
}

double OfflatticeMonteCarlo::concentration() const
{
    return nOfflatticeParticles()/solver().freeVolume();
}

bool OfflatticeMonteCarlo::hasDiscreteParticles() const
{
    return true;
}

uint OfflatticeMonteCarlo::numberOfParticles() const
{
    return nOfflatticeParticles();
}

void OfflatticeMonteCarlo::insertRandomParticle()
{
    const double zMin = solver().heights().min();
    const double &h = solver().confiningSurfaceEvent().height();

    double x0;
    double y0;
    double z0;

    do
    {
        double x0Raw = rng.uniform()*solver().length();
        double y0Raw = rng.uniform()*solver().width();
        z0 = zMin + rng.uniform()*(h - zMin);

        //            //it is rarely the case that the boundary transformations
        //            //touch a point inside the solver cube.
        //            x0 = solver().boundaryTransform(x0Raw, y0Raw, z0, 0);
        //            y0 = solver().boundaryTransform(x0, y0Raw, z0, 1);

        solver().boundaryContinousTransform(x0, y0, x0Raw, y0Raw, z0);

    } while(solver().isBlockedPosition(x0, y0, z0));

    insertParticle(x0, y0, z0);

}

void OfflatticeMonteCarlo::removeRandomParticle()
{
    const uint n = rng.uniform()*numberOfParticles();

    removeParticle(n);
}

void OfflatticeMonteCarlo::diffuse(const double dt)
{
    if (dt < 0 || fabs(dt) < 1E-15)
    {
        return;
    }

    double x1, y1, z1;
    double x0, y0, z0;

    //we use DUnscaled because dt is in unscaled units
    const double D = DUnscaled();
    const double prefac = sqrt(2*D*dt);

    for (uint n = 0; n < nOfflatticeParticles(); ++n)
    {
        x0 = m_particlePositions(0, n);
        y0 = m_particlePositions(1, n);
        z0 = m_particlePositions(2, n);

        BADAssBool(!solver().isBlockedPosition(x0, y0, z0));

        //        m_F(0, n) = 0;
        //        m_F(1, n) = 0;
        //        m_F(2, n) = solver().confiningSurfaceEvent().diffusionDrift(x0, y0, z0);


        do
        {
            double dx = prefac*rng.normal();// + D*m_F(0, n)*dt;
            double dy = prefac*rng.normal();// + D*m_F(1, n)*dt;
            double dz = prefac*rng.normal();// + D*m_F(2, n)*dt;

            //            x1 = solver().boundaryTransform(x0, y0, z0, dx, 0);
            //            y1 = solver().boundaryTransform(x0, y0, z0, dy, 1);

            solver().boundaryContinousTransform(x1, y1, x0 + dx, y0 + dy, z0);
            z1 = z0 + dz;

            if (solver().surfaceDim() == 1)
            {
                y1 = y0;
            }



        } while(solver().isBlockedPosition(x1, y1, z1));

        if (solver().confiningSurfaceEvent().acceptDiffusionMove(x0, y0, z0, x1, y1, x1))
        {
            m_particlePositions(0, n) = x1;
            m_particlePositions(1, n) = y1;
            m_particlePositions(2, n) = z1;

            m_accepted++;
        }

        m_trials++;
    }
}

void OfflatticeMonteCarlo::diffuseFull(const double dtFull)
{
    const uint N = dtFull/maxdt();

    //    cout << setprecision(16) << fixed << D() << " " << dtFull << " " << sqrt(2*D()*dtFull) << endl;

    for (uint i = 0; i < N; ++i)
    {
        diffuse(maxdt());
    }

    diffuse(dtFull - N*maxdt());

    calculateLocalRates();
}

void OfflatticeMonteCarlo::removeParticle(const uint n)
{
    //    m_particlePositions.shed_col(n);
    //    m_F.shed_col(n);

    //    m_localRates.shed_slice(n);

    m_nParticles--;

    for (uint dim = 0; dim < 3; ++dim)
    {
        m_particlePositions(dim, n) = m_particlePositions(dim, m_nParticles);
        m_F(dim, n) = m_F(dim, m_nParticles);
    }

    m_localRates.slice(n) = m_localRates.slice(m_nParticles);

}

void OfflatticeMonteCarlo::insertParticle(const double x, const double y, const double z)
{

    if (m_particlePositions.n_cols <= m_nParticles)
    {
        m_particlePositions.resize(3, m_nParticles*2);

        m_F.resize(3, m_nParticles*2);

        m_localRates.resize(solver().length(), solver().width(), m_nParticles*2);

        m_localRatesForSite.resize(m_nParticles*2);
    }

    m_particlePositions(0, m_nParticles) = x;
    m_particlePositions(1, m_nParticles) = y;
    m_particlePositions(2, m_nParticles) = z;

    m_nParticles++;

    //    for (uint x = 0; x < solver().length(); ++x)
    //    {
    //        for (uint y = 0; y < solver().width(); ++y)
    //        {
    //            m_localRates(x, y, nOfflatticeParticles() - 1) = calculateLocalRateOverD(x, y, nOfflatticeParticles() - 1);
    //        }
    //    }

}

void OfflatticeMonteCarlo::initializeParticleMatrices(const uint N)
{
    while (m_nParticles < N)
    {
        insertRandomParticle();
    }
}

void OfflatticeMonteCarlo::scan(const uint n, const uint dim, const double dr, const uint maxSteps)
{
    const double &x0 = particlePositions(0, n);
    const double &y0 = particlePositions(1, n);
    const double &z0 = particlePositions(2, n);

    uint c = 0;

    while (solver().isBlockedPosition(x0, y0, z0) && c < maxSteps)
    {
        m_particlePositions(dim, n) = solver().boundaryContinousTransformSingle(particlePositions(0, n),
                                                                                particlePositions(1, n),
                                                                                particlePositions(2, n),
                                                                                dim, dr);

        c++;
    }
}

void OfflatticeMonteCarlo::scanForDisplacement(const uint n, uint &dim, double &delta, const double stepSize)
{
    m_scanOriginalPositions(0) = m_particlePositions(0, n);
    m_scanOriginalPositions(1) = m_particlePositions(1, n);
    m_scanOriginalPositions(2) = m_particlePositions(2, n);

    uint c = 0;
    for (uint dim = 0; dim < 3; ++dim)
    {
        for (int direction = -1; direction <= 1; direction += 2)
        {
            scan(n, dim, direction*stepSize);
            m_scanDeltas(c) = particlePositions(dim, n) - m_scanOriginalPositions(dim);
            m_scanAbsDeltas(c) = fabs(m_scanDeltas(c));
            m_particlePositions(dim, n) = m_scanOriginalPositions(dim);
            c++;
        }
    }

    uint minLoc;

    m_scanAbsDeltas.min(minLoc);
    delta = m_scanDeltas(minLoc);
    dim = minLoc/2;

}

void OfflatticeMonteCarlo::dumpDiffusingParticles(const uint frameNumber, const string path) const
{
    const double &h = solver().confiningSurfaceEvent().height();
    const int zMin = solver().heights().min();

    lammpswriter writer(5, "cavitydiff", path);
    writer.setSystemSize(solver().length(), solver().width(), h, 0, 0, zMin);
    writer.initializeNewFile(frameNumber);

    if (nOfflatticeParticles() == 0)
    {
        writer << 0 << 0 << 0 << zMin << 0;
    }

    for (uint n = 0; n < nOfflatticeParticles(); ++n)
    {
        const double &x = m_particlePositions(0, n);
        const double &y = m_particlePositions(1, n);
        const double &z = m_particlePositions(2, n);

        BADAssBool(!solver().isBlockedPosition(x, y, z));

        writer << 0
               << x
               << y
               << z
               << totalParticleDepositionRate(n);
    }

    writer.finalize();
}

void OfflatticeMonteCarlo::clearDiffusingParticles()
{
    //    m_particlePositions.clear();
    //    m_F.clear();
    //    m_localRates.clear();

    m_nParticles = 0;
}

void OfflatticeMonteCarlo::calculateLocalRates()
{
    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            solver().surfaceReaction(x, y).setDepositionRate(0);
        }
    }

    const int &l = m_depositionBoxHalfSize;

    int xTrans;
    int yTrans;

    DissolutionDeposition *r;
    double r2;
    double localRate;

    const int hmax = solver().heights().max();

    const double D = DScaled();

    //    lammpswriter writer(4, "astar", "/tmp");
    //    lammpswriter surf(3, "surf", "/tmp");
    //    surf.setSystemSize(2*l+1, 2*l+1, 2*l+1);
    //    writer.setSystemSize(2*l+1, 2*l+1, 2*l+1);


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
                        m_world->MarkPosition(xw, yw, zwIter, true);
                        //                        surf << xw << yw << zwIter;
                    }

                    m_pathFindingJazzes[m_nPathFinds]->xTrans = xTrans;
                    m_pathFindingJazzes[m_nPathFinds]->yTrans = yTrans;
                    m_pathFindingJazzes[m_nPathFinds]->xEnd = xw;
                    m_pathFindingJazzes[m_nPathFinds]->yEnd = yw;
                    m_pathFindingJazzes[m_nPathFinds]->zEnd = zw;
                    m_nPathFinds++;

                    //                    r2 = dz2 + (xscan + dx)*(xscan + dx) + (yscan + dy)*(yscan+dy);

                    //                    BADAssClose(r2, solver().closestSquareDistance(xTrans, yTrans, h+1, xp, yp, zp), 1E-3);

                    //                    r = &solver().surfaceReaction(xTrans, yTrans);
                    //                    localRate = calculateLocalRateOverD(r2);
                    //                    m_localRates(xTrans, yTrans, n) = localRate;
                    //                    r->setDepositionRate(r->depositionRate() + localRate*D);
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

            r2 = dx*dx + dy*dy + dz*dz;
            const SearchNode *crumb = m_pathFinder->findPath(l, l, l, xw, yw, zw);
            //            writer << 3+i << crumb->position->X << crumb->position->Y << crumb->position->Z;
            if (crumb == nullptr)
            {
                continue;
            }

            while (crumb->next != nullptr)
            {
                r2 += crumb->GetDistanceSquared(crumb->next->x, crumb->next->y, crumb->next->z);
                //                writer << 3+i << crumb->position->X << crumb->position->Y << crumb->position->Z;

                crumb = crumb->next;
            }
            //            writer << 3+i << crumb->position->X << crumb->position->Y << crumb->position->Z;

            r = &solver().surfaceReaction(xTrans, yTrans);
            localRate = calculateLocalRateOverD(r2);
            m_localRates(xTrans, yTrans, n) = localRate;
            r->setDepositionRate(r->depositionRate() + localRate*D);
        }



        //        if (surf.valueCounter() == 0)
        //        {
        //            surf << 0 << 0 << 0;
        //        }

        //        writer.finalize();
        //        surf.finalize();
    }

    return;

    //    DissolutionDeposition *r;
    //    for (uint _x = 0; _x < solver().length(); ++_x)
    //    {
    //        for (uint _y = 0; _y < solver().width(); ++_y)
    //        {
    //            r = &solver().surfaceReaction(_x, _y);

    //            const double newDepositionRate = depositionRate(_x, _y);

    //            r->setDepositionRate(newDepositionRate);
    ////            r->changeRate(newDepositionRate + r->dissolutionRate());
    //        }
    //    }

    //    for (uint x = 0; x < solver().length(); ++x)
    //    {
    //        for (uint y = 0; y < solver().width(); ++y)
    //        {
    //            for (uint n = 0; n < nOfflatticeParticles(); ++n)
    //            {
    //                m_localRates(x, y, n) = calculateLocalRateOverD(x, y, n);
    //            }
    //        }
    //    }

    //    selectDepositionReactants();
}

double OfflatticeMonteCarlo::totalParticleDepositionRate(const uint n) const
{
    double r = 0;

    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            r += m_localRates(x, y, n);
        }
    }

    return r;
}

void OfflatticeMonteCarlo::selectDepositionReactants()
{
    uint xSelected, ySelected;
    for (uint n = 0; n < nOfflatticeParticles(); ++n)
    {
        selectDepositionReactant(xSelected, ySelected, n);

        m_particleDepositionLocations(n, 0) = xSelected;
        m_particleDepositionLocations(n, 1) = ySelected;
    }
}

void OfflatticeMonteCarlo::selectDepositionReactant(uint &xSelected, uint &ySelected, const uint n)
{
    vec accu(solver().length()*solver().width());

    double totalRate = 0;
    uint c = 0;
    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            totalRate += localRates(x, y, n);

            accu(c) = totalRate;

            c++;
        }
    }

    const uint choice = chooseFromTotalRate(accu, totalRate);

    xSelected = choice / solver().width();
    ySelected = choice % solver().width();
}

bool OfflatticeMonteCarlo::isInLineOfSight(const uint n, const uint x, const uint y) const
{
    //slow as ****

    const int z = solver().height(x, y) + 1;

    double xp = particlePositions(n, 0);
    double yp = particlePositions(n, 1);
    double zp = particlePositions(n, 2);

    const double dx = x - xp;
    const double dy = y - yp;
    const double dz = z - zp;

    const double r = sqrt(dx*dx + dy*dy + dz*dz);

    const double rxy = sqrt(dx*dx + dy*dy);

    const double xyAngle = atan2(dy, dx);
    const double yzAngle = atan2(dz, rxy);

    const double _dr = 0.01;
    const double _dz = _dr*sin(yzAngle);

    const double _drxy = _dz/dz*rxy;

    const double _dx = _drxy*cos(xyAngle);
    const double _dy = _drxy*sin(xyAngle);

    const uint nTraceSteps = r/_dr;

    for (uint nTrace = 0; nTrace < nTraceSteps; ++nTrace)
    {
        xp += _dx;
        yp += _dy;
        zp += _dz;

        if (solver().isBlockedPosition(xp, yp, zp))
        {
            return false;
        }
    }

    return true;

}

void OfflatticeMonteCarlo::initialize()
{
    m_accepted = 0;
    m_trials = 0;
}

void OfflatticeMonteCarlo::reset()
{
    const double &dtFull = solver().currentTimeStep();

    diffuseFull(dtFull);
}

void OfflatticeMonteCarlo::initializeObserver(const Subjects &subject)
{
    (void) subject;

    //already initialized
    if (numberOfParticles() == 0)
    {
        const double V = solver().freeVolume();

        const uint N = V*solver().concentration() + rng.uniform();

        initializeParticleMatrices(N);
    }

    calculateLocalRates();
}

void OfflatticeMonteCarlo::execute()
{
}
