#include "offlatticemontecarlo.h"

#include "../confiningsurface/confiningsurface.h"
#include "fluxboundaryreaction.h"

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


void OfflatticeMonteCarlo::dump(const uint frameNumber, const string path, const string ext) const
{
    Diffusion::dump(frameNumber, path, ext);
    dumpDiffusingParticles(frameNumber, path, ext);
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
                if (!solver().fluxBoundaryDeposition())
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
            }

            //Add a particle on a unit sphere around the site;
            else
            {
                int dx, dy, dz;
                int x1, y1, z1;

                const int z = solver().height(x, y) + 1;

                solver().getRandomSolutionSite(x, y, z, dx, dy, dz);

                solver().boundaryLatticeTransform(x1, y1, x + dx, y + dy, z);

                if (!solver().isOutsideBoxContinuous(x1, y1))
                {
                    z1 = z + dz;

                    insertParticle(x1, y1, z1);

                    m_lockedParticles.push_back(nOfflatticeParticles() - 1);
                }
            }
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

                BADAssBool(!solver().isOutsideBoxContinuous(x0, y0));

                scanForDisplacement(n, dim, delta);
                m_particlePositions(dim, n) += delta;

                BADAssBool(!solver().isBlockedPosition(x0, y0, z0));
                BADAssBool(!solver().isOutsideBoxContinuous(x0, y0));
            }
        }
    }
    
    else
    {
        const double &hc = solver().confiningSurfaceEvent().height();

        vector<uint> lostParticles;
        for (uint n = 0; n < nOfflatticeParticles(); ++n)
        {
            const double x0 = particlePositions(0, n);
            const double y0 = particlePositions(1, n);
            const double z0 = particlePositions(2, n);

            //this means that the particle has been locked inside
            //the surface when it moved
            if (z0 > hc - 1)
            {
                const int &hUnder = solver().height(round(x0), round(y0));

                //if there is no room for it directly below..
                if (hc - hUnder <= 2)
                {
                    //..we remove the particle and insert it at a
                    //random spot
                    lostParticles.push_back(n);
                }

                else if (hUnder == hc - 2)
                {
                    m_particlePositions(2, n) = hUnder + 1;
                }

                else
                {
                    //..else we shift it down to be just below the
                    //confining surface
                    m_particlePositions(2, n) = hc - 1.0;
                }
            }
        }

        //remove in decreasing order to avoid messing things up
        for (int i = int(lostParticles.size() - 1); i >= 0; --i)
        {
            removeParticle(lostParticles.at(i));
        }

        insertRandomParticles(lostParticles.size());

#ifndef NDEBUG
        for (uint n = 0; n < nOfflatticeParticles(); ++n)
        {
            BADAss(particlePositions(2, n), <=, hc - 1,
                   "illegal translation", [&] ()
            {
                const double x = particlePositions(0, n);
                const double y = particlePositions(1, n);
                const double z = particlePositions(2, n);
                const uint X = round(x);
                const uint Y = round(y);
                const int &h = solver().height(X, Y);
                BADAssSimpleDump(hc, x, y, z, X, Y, h);
            });
        }
#endif
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

void OfflatticeMonteCarlo::executeFluxBoundaryReaction(const uint x, const uint y, const double z)
{
    insertParticle(x, y, z);
    m_lockedParticles.push_back(m_nParticles - 1);
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
    return nOfflatticeParticles()/volume();
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
    double x0;
    double y0;
    double z0;

    findRandomPosition(x0, y0, z0);
    insertParticle(x0, y0, z0);
}

void OfflatticeMonteCarlo::removeRandomParticle()
{
    const uint n = rng.uniform()*numberOfParticles();

    removeParticle(n);
}

void OfflatticeMonteCarlo::diffuse(const double dt, const double prefac)
{
    if (nOfflatticeParticles() == 0 || dt < 1E-15)
    {
        return;
    }

    for (uint n = 0; n < nOfflatticeParticles(); ++n)
    {
        if (std::find(m_lockedParticles.begin(), m_lockedParticles.end(), n) != m_lockedParticles.end())
        {
            continue;
        }

        if (std::find(m_removalQueue.begin(), m_removalQueue.end(), n) != m_removalQueue.end())
        {
            continue;
        }

        diffuseSingleParticle(n, dt, prefac);
    }
}

void OfflatticeMonteCarlo::diffuseFull(const double dtFull)
{
    const uint N = dtFull/maxdt();

    const double maxPrefac = sqrt(2*DUnscaled()*maxdt());
    for (uint i = 0; i < N; ++i)
    {
        diffuse(maxdt(), maxPrefac);
    }

    diffuse(dtFull - N*maxdt());
}

void OfflatticeMonteCarlo::diffuseSingleParticle(const uint n, const double dt, const double prefac)
{
    (void) dt;
    double x1, y1, z1;

    const double &x0 = m_particlePositions(0, n);
    const double &y0 = m_particlePositions(1, n);
    const double &z0 = m_particlePositions(2, n);

    BADAssBool(!solver().isBlockedPosition(x0, y0, z0), "hmm", [&] ()
    {
        dump(cycle() + 1);
        int h = solver().height(round(x0), round(y0));
        double hc = solver().confiningSurfaceEvent().height();
        BADAssSimpleDump(cycle(), x0, y0, z0, h, hc);
    });

    BADAssBool(!solver().isOutsideBoxContinuous(x0, y0));

    //        m_F(0, n) = 0;
    //        m_F(1, n) = 0;
    //        m_F(2, n) = solver().confiningSurfaceEvent().diffusionDrift(x0, y0, z0);

    //if particle is squeezed we would wait
    //forever for a dz=0 trial position.
    bool dzZero = false;
    const bool eps = 1E-3;
    if (fabs(solver().confiningSurfaceEvent().height() - z0) < 1 + eps &&
            fabs(solver().height(round(x0), round(y0)) - z0) < 1 + eps)
    {
        dzZero = true;
    }

    double dx;
    double dy;
    double dz;

    uint N = 0;
    const uint NMax = 1000000;

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

        if (solver().isOutsideBoxContinuous(x1, y1))
        {
            m_removalQueue.push_back(n);
            return;
        }

        N++;

        if (N > NMax)
        {
            cout << cycle() + 1 << " " << n << " " << x0 << " " << y0 << " " << z0 << endl;
            dump(cycle() + 1);

            throw std::runtime_error("No path available for diffusing particle.");
        }

    } while(solver().isBlockedPosition(x1, y1, z1));

    if (solver().confiningSurfaceEvent().acceptDiffusionMove(x0, y0, z0, x1, y1, x1))
    {
        m_particlePositions(0, n) = x1;
        m_particlePositions(1, n) = y1;
        m_particlePositions(2, n) = z1;
    }
}

void OfflatticeMonteCarlo::removeParticle(const uint n)
{
    m_nParticles--;

    BADAss(n, <=, m_nParticles);

    if (n != m_nParticles)
    {

        for (uint dim = 0; dim < 3; ++dim)
        {
            m_particlePositions(dim, n) = m_particlePositions(dim, m_nParticles);
            m_F(dim, n) = m_F(dim, m_nParticles);
        }

        m_localRates.slice(n) = m_localRates.slice(m_nParticles);
    }

}

void OfflatticeMonteCarlo::insertParticle(const double x, const double y, const double z)
{
    BADAssBool(!solver().isBlockedPosition(x, y, z));
    BADAssBool(!solver().isOutsideBoxContinuous(x, y));

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

        if (solver().isOutsideBoxContinuous(x0, y0))
        {
            break;
        }
    }
}

void OfflatticeMonteCarlo::scanForDisplacement(const uint n, uint &dim, double &delta, const double stepSize)
{
    m_scanOriginalPositions(0) = m_particlePositions(0, n);
    m_scanOriginalPositions(1) = m_particlePositions(1, n);
    m_scanOriginalPositions(2) = m_particlePositions(2, n);

    bool noWhereToGo = true;

    uint c = 0;
    for (uint dim = 0; dim < 3; ++dim)
    {
        for (int direction = -1; direction <= 1; direction += 2)
        {
            scan(n, dim, direction*stepSize);

            bool available = true;

            if (solver().isOutsideBoxContinuous(particlePositions(0, n),
                                                particlePositions(1, n)))
            {
                available = false;
            }

            else
            {
                if (solver().isBlockedPosition(particlePositions(0, n),
                                               particlePositions(1, n),
                                               particlePositions(2, n)))
                {
                    available = false;
                }
            }

            if (available)
            {
                const double delta = particlePositions(dim, n) - m_scanOriginalPositions(dim);
                noWhereToGo = false;
                m_scanDeltas(c) = delta;
            }

            else
            {
                m_scanDeltas(c) = numeric_limits<double>::max();
            }

            m_scanAbsDeltas(c) = fabs(m_scanDeltas(c));
            m_particlePositions(dim, n) = m_scanOriginalPositions(dim);
            c++;
        }
    }

    if (noWhereToGo)
    {
//        if (noAvailablePositions())
//        {
//            return;
//        }

        double x0;
        double y0;
        double z0;

        findRandomPosition(x0, y0, z0);

        BADAssBool(!solver().isBlockedPosition(x0, y0, z0));
        BADAssBool(!solver().isOutsideBoxContinuous(x0, y0));

        m_particlePositions(0, n) = x0;
        m_particlePositions(1, n) = y0;
        m_particlePositions(2, n) = z0;

        delta = 0;
        dim = 0;
        return;
    }

    uint minLoc = 0;

    m_scanAbsDeltas.min(minLoc);
    delta = m_scanDeltas(minLoc);
    dim = minLoc/2;
}

bool OfflatticeMonteCarlo::noAvailablePositions() const
{
    const double &h = solver().confiningSurfaceEvent().height();

    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            if (h - solver().height(x, y) - 2 > 0)
            {
                return false;
            }
        }
    }

    return true;
}

void OfflatticeMonteCarlo::dumpDiffusingParticles(const uint frameNumber, const string path, const string ext) const
{
    const double &h = solver().confiningSurfaceEvent().height();
    const int zMin = solver().heights().min();

    lammpswriter writer(5, "cavitydiff" + ext, path);
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

        if (solver().isBlockedPosition(xp, yp, zp) || solver().isOutsideBoxContinuous(xp, yp))
        {
            return false;
        }
    }

    return true;

}

void OfflatticeMonteCarlo::releaseLockedParticles()
{
    m_lockedParticles.clear();
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

void OfflatticeMonteCarlo::findRandomPosition(double &x0, double &y0, double &z0) const
{
    const double zMin = solver().heights().min() + 1;
    const double &h = solver().confiningSurfaceEvent().height() - 1;

    const uint nMax = 10000;
    uint n = 0;

    do
    {
        double x0Raw = -0.5 + rng.uniform()*solver().length();
        double y0Raw = -0.5 + rng.uniform()*solver().width();
        z0 = zMin + rng.uniform()*(h - zMin);

        solver().boundaryContinousTransform(x0, y0, x0Raw, y0Raw, z0);

        n++;

        if (n > nMax)
        {
            throw std::runtime_error("random position asked for, but no positions are legal.");
        }

    } while(solver().isBlockedPosition(x0, y0, z0));
}

double OfflatticeMonteCarlo::volume() const
{
    return solver().volume();
}

void OfflatticeMonteCarlo::initialize()
{
    releaseLockedParticles();
    m_dumpCounter = 0;
}

void OfflatticeMonteCarlo::reset()
{
    const double &dtFull = solver().currentTimeStep();

    m_removalQueue.clear();

    diffuseFull(dtFull);

    releaseLockedParticles();

    std::sort(m_removalQueue.begin(), m_removalQueue.end());
    for (int i = m_removalQueue.size() - 1; i >= 0; --i)
    {
        const uint n = m_removalQueue.at(i);
        removeParticle(n);
    }

    calculateLocalRatesAndUpdateDepositionRates();
}

void OfflatticeMonteCarlo::initializeObserver(const Subjects &subject)
{
    (void) subject;

    //already initialized
    if (numberOfParticles() == 0)
    {
        const double V = volume();

        BADAss(V, >=, 0);

        const uint N = V*solver().concentration() + rng.uniform();

        insertRandomParticles(N);
    }

    calculateLocalRatesAndUpdateDepositionRates();
}
