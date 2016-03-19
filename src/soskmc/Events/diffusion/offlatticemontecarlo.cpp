#include "offlatticemontecarlo.h"

#include "../confiningsurface/confiningsurface.h"
#include "concentrationboundaryreaction.h"

OfflatticeMonteCarlo::OfflatticeMonteCarlo(SOSSolver &solver,
                                           const double maxdt,
                                           string type,
                                           string unit,
                                           bool hasOutput,
                                           bool storeValue) :
    Diffusion(solver, type, unit, hasOutput, storeValue),
    m_localRates(solver.length(), solver.width(), 100),
    m_maxdt(maxdt),
    m_particlePositions(3, 100),
    m_F(3, 100, fill::zeros),
    m_accuRatesForSite(100),
    m_nParticles(0)
{

}

OfflatticeMonteCarlo::~OfflatticeMonteCarlo()
{

}


void OfflatticeMonteCarlo::dump(const uint frameNumber, const string path) const
{
    Diffusion::dump(frameNumber, path);
    dumpDiffusingParticles(frameNumber, path);
}

bool OfflatticeMonteCarlo::countPaths() const
{
    return true;
}

void OfflatticeMonteCarlo::notifyObserver(const Subjects &subject)
{
    if (!hasStarted())
    {
        return;
    }

    if (subject == Subjects::SOLVER)
    {
        if (solver().currentSurfaceChange().type == ChangeTypes::Single)
        {
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
                    Rtot += localRates(x, y, n);
                    m_accuRatesForSite(n) = Rtot;
                }

                const uint N = chooseFromTotalRate(m_accuRatesForSite.memptr(), nOfflatticeParticles(), Rtot);

                BADAss(localRates(x, y, N), !=, 0);

                removeParticle(N);
            }

            //Add a particle on a unit sphere around the site;
            else
            {
                int dx, dy, dz;
                int x1, y1, z1;

                const int &z = solver().height(x, y) + 1;

                solver().getRandomSolutionSite(x, y, z, dx, dy, dz);

                solver().boundaryLatticeTransform(x1, y1, x + dx, y + dy, z);

                z1 = z + dz;

                insertParticle(x1, y1, z1);

                m_particleIsLocked = true;
                m_lockedParticle = nOfflatticeParticles() - 1;
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
                }
            }
        }
    }
    
    else
    {
        const double &hc = solver().confiningSurfaceEvent().height();

        vector<uint> lostParticles;
        for (uint n = 0; n < nOfflatticeParticles(); ++n)
        {
            const double &x0 = particlePositions(0, n);
            const double &y0 = particlePositions(1, n);
            const double &z0 = particlePositions(2, n);

            //this means that the particle has been locked inside
            //the surface when it moved
            if (z0 > hc - 1)
            {
                //if there is no room for it directly below..
                if (solver().height(round(x0), round(y0)) == hc - 1)
                {
                    //..we remove the particle and insert it at a
                    //random spot
                    lostParticles.push_back(n);
                }
                else
                {
                    //..else we shift it down to be just below the
                    //confining surface
                    m_particlePositions(2, n) = hc - 1;
                }
            }
        }

        //remove in decreasing order to avoid messing things up
        for (int i = int(lostParticles.size() - 1); i >= 0; ++i)
        {
            removeParticle(lostParticles.at(i));
        }

        insertRandomParticles(lostParticles.size());
    }
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
    if (nOfflatticeParticles() == 0 || dt < 1E-15)
    {
        return;
    }

    //we use DUnscaled because dt is in unscaled units
    const double prefac = sqrt(2*DUnscaled()*dt);

    for (uint n = 0; n < nOfflatticeParticles(); ++n)
    {
        if (m_particleIsLocked)
        {
            if (n == m_lockedParticle)
            {
                continue;
            }

        }

        diffuseSingleParticle(n, dt, prefac);
    }
}

void OfflatticeMonteCarlo::diffuseFull(const double dtFull)
{
    const uint N = dtFull/maxdt();

    for (uint i = 0; i < N; ++i)
    {
        diffuse(maxdt());
    }

    diffuse(dtFull - N*maxdt());

    calculateLocalRatesAndUpdateDepositionRates();
}

void OfflatticeMonteCarlo::diffuseSingleParticle(const uint n, const double dt, const double prefac)
{
    (void) dt;
    double x1, y1, z1;

    const double &x0 = m_particlePositions(0, n);
    const double &y0 = m_particlePositions(1, n);
    const double &z0 = m_particlePositions(2, n);

    BADAssBool(!solver().isBlockedPosition(x0, y0, z0));

    //        m_F(0, n) = 0;
    //        m_F(1, n) = 0;
    //        m_F(2, n) = solver().confiningSurfaceEvent().diffusionDrift(x0, y0, z0);

    //if particle is squeezed we would wait
    //forever for a dz=0 trial position.
    bool dzZero = false;
    if (solver().confiningSurfaceEvent().height() == z0 + 1 &&
            solver().height(round(x0), round(y0)) == z0 - 1)
    {
        dzZero = true;
    }

    double dx;
    double dy;
    double dz;

    do
    {
        dx = prefac*rng.normal();// + D*m_F(0, n)*dt;
        dy = prefac*rng.normal();// + D*m_F(1, n)*dt;

        solver().boundaryContinousTransform(x1, y1, x0 + dx, y0 + dy, z0);

        if (!dzZero)
        {
            dz = prefac*rng.normal();// + D*m_F(2, n)*dt;
            z1 = z0 + dz;
        }

        else
        {
            z1 = z0;
        }

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

void OfflatticeMonteCarlo::removeParticle(const uint n)
{
    m_nParticles--;

    for (uint dim = 0; dim < 3; ++dim)
    {
        m_particlePositions(dim, n) = m_particlePositions(dim, m_nParticles);
        m_F(dim, n) = m_F(dim, m_nParticles);
    }

    BADAss(n, <=, m_nParticles);

    if (n != m_nParticles)
    {
        m_localRates.slice(n) = m_localRates.slice(m_nParticles);
    }

}

void OfflatticeMonteCarlo::insertParticle(const double x, const double y, const double z)
{

    if (m_particlePositions.n_cols <= m_nParticles)
    {
        m_particlePositions.resize(3, m_nParticles*2);

        m_F.resize(3, m_nParticles*2);

        m_localRates.resize(solver().length(), solver().width(), m_nParticles*2);

        m_accuRatesForSite.resize(m_nParticles*2);
    }

    m_particlePositions(0, m_nParticles) = x;
    m_particlePositions(1, m_nParticles) = y;
    m_particlePositions(2, m_nParticles) = z;

    m_nParticles++;

}

void OfflatticeMonteCarlo::insertRandomParticles(const uint N)
{
    uint nParticles = m_nParticles + N;
    while (m_nParticles < nParticles)
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

    uint minLoc = 0;

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
    m_nParticles = 0;
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

void OfflatticeMonteCarlo::releaseLockedParticle()
{
    m_particleIsLocked = false;
}

void OfflatticeMonteCarlo::resetLocalRates(const uint n)
{
    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            m_localRates(x, y, n) = 0;
        }
    }
}

void OfflatticeMonteCarlo::initialize()
{
    m_particleIsLocked = false;

    m_accepted = 0;
    m_trials = 0;
}

void OfflatticeMonteCarlo::reset()
{
    const double &dtFull = solver().currentTimeStep();

    diffuseFull(dtFull);

    releaseLockedParticle();
}

void OfflatticeMonteCarlo::initializeObserver(const Subjects &subject)
{
    (void) subject;

    //already initialized
    if (numberOfParticles() == 0)
    {
        const double V = solver().freeVolume();

        BADAss(V, >=, 0);

        const uint N = V*solver().concentration() + rng.uniform();

        insertRandomParticles(N);
    }

    calculateLocalRatesAndUpdateDepositionRates();
}
