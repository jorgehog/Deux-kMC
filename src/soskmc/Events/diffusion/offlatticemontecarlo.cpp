#include "offlatticemontecarlo.h"

#include "../confiningsurface/confiningsurface.h"
#include "concentrationboundaryreaction.h"

#include "dissolutiondeposition.h"

OfflatticeMonteCarlo::OfflatticeMonteCarlo(SOSSolver &solver,
                                           const double maxdt,
                                           string type,
                                           string unit,
                                           bool hasOutput,
                                           bool storeValue) :
    Diffusion(solver, type, unit, hasOutput, storeValue),
    m_mutexSolver(solver),
    m_maxdt(maxdt)
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

void OfflatticeMonteCarlo::registerHeightChange(const uint x,
                                                const uint y,
                                                const int value,
                                                std::vector<DissolutionDeposition *> &affectedSurfaceReactions,
                                                const uint n)
{
    (void) affectedSurfaceReactions;
    (void) n;

    if (!hasStarted())
    {
        return;
    }

    //Remove particle based on the probability of it being the deposited
    if (value == 1)
    {
        vector<double> localRatesForSite(nOfflatticeParticles());

        double Rtot = 0;
        for (uint n = 0; n < nOfflatticeParticles(); ++n)
        {
            //use old rates to calculate probability of depositing
            localRatesForSite.at(n) = localRates(x, y, n);
            Rtot += localRatesForSite.at(n);
        }

        double R = rng.uniform()*Rtot;

        uint N = binarySearchForInterval(R, localRatesForSite);

        removeParticle(N);
    }

    //Add a particle on a unit sphere around the site;
    else
    {
        double dx, dy, dz;
        double x1, y1, z1;

        const int &z = solver().height(x, y) + 1;

        do
        {
            double theta = rng.uniform()*datum::pi;
            double phi = rng.uniform()*datum::pi*2;

            dx = sin(theta)*cos(phi);
            dy = sin(theta)*sin(phi);
            dz = cos(theta);

            x1 = solver().boundaryTransform(x, y, z, dx, 0);
            y1 = solver().boundaryTransform(x, y, z, dy, 1);
            z1 = z + dz;

            BADAssClose(1, dx*dx + dy*dy + dz*dz, 1E-3);


        } while (solver().isBlockedPosition(x1, y1, z1));

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
            particlePositions(dim, n) += delta;

            for (uint _x = 0; _x < solver().length(); ++_x)
            {
                for (uint _y = 0; _y < solver().width(); ++_y)
                {
                    m_localRates(_x, _y, n) = calculateLocalRate(_x, _y, n);
                }
            }
        }

        else
        {
            //else we just recalulate the rate for the new height
            m_localRates(x, y, n) = calculateLocalRate(x, y, n);
        }
    }

#ifndef NDEBUG
    for (uint _x = 0; _x < solver().length(); ++_x)
    {
        for (uint _y = 0; _y < solver().width(); ++_y)
        {
            for (uint n = 0; n < nOfflatticeParticles(); ++n)
            {
                BADAssClose(m_localRates(_x, _y, n), calculateLocalRate(_x, _y, n), 1E-3, "what", [&] ()
                {
                    BADAssSimpleDump(_x, _y, n);
                });
            }
        }
    }
#endif

    DissolutionDeposition *r;
    for (uint _x = 0; _x < solver().length(); ++_x)
    {
        for (uint _y = 0; _y < solver().width(); ++_y)
        {
            r = &m_mutexSolver.surfaceReaction(_x, _y);

            const double newDepositionRate = depositionRate(_x, _y);

            r->setDepositionRate(newDepositionRate);
            r->changeRate(newDepositionRate + r->dissolutionRate());
        }
    }

}

double OfflatticeMonteCarlo::depositionRate(const uint x, const uint y) const
{

    //if there is no site available between the confining wall and the surface,
    //the rate is zero.
    if (solver().isBlockedPosition(x, y, solver().height(x, y) + 1))
    {
        return 0.0;
    }

    double R = 0;

    for (uint n = 0; n < nOfflatticeParticles(); ++n)
    {
        const double Rn = localRates(x, y, n);

        BADAssClose(Rn, calculateLocalRate(x, y, n), 1E-3);

        R += Rn;
    }

    return R/solver().concentration();
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

void OfflatticeMonteCarlo::diffuse(const double dt)
{
    if (dt < 0 || fabs(dt) < 1E-15)
    {
        return;
    }

    double x1, y1, z1;
    double x0, y0, z0;

    for (uint n = 0; n < nOfflatticeParticles(); ++n)
    {
        x0 = m_particlePositions(0, n);
        y0 = m_particlePositions(1, n);
        z0 = m_particlePositions(2, n);

        BADAssBool(!solver().isBlockedPosition(x0, y0, z0));

        m_F(2, n) = solver().confiningSurfaceEvent().diffusionDrift(x0, y0, z0);

        do
        {
            double dx = sqrt(2*D()*dt)*rng.normal() + D()*m_F(0, n)*dt;
            double dy = sqrt(2*D()*dt)*rng.normal() + D()*m_F(1, n)*dt;
            double dz = sqrt(2*D()*dt)*rng.normal() + D()*m_F(2, n)*dt;

            x1 = solver().boundaryTransform(x0, y0, z0, dx, 0);
            y1 = solver().boundaryTransform(x0, y0, z0, dy, 1);
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

    for (uint i = 0; i < N; ++i)
    {
        diffuse(maxdt());
    }

    diffuse(dtFull - N*maxdt());

    calculateLocalRates();
}

void OfflatticeMonteCarlo::removeParticle(const uint n)
{
    m_particlePositions.shed_col(n);
    m_F.shed_col(n);

    m_localRates.shed_slice(n);
}

void OfflatticeMonteCarlo::insertParticle(const double x, const double y, const double z)
{
    m_particlePositions.resize(3, nOfflatticeParticles() + 1);

    m_particlePositions(0, nOfflatticeParticles() - 1) = x;
    m_particlePositions(1, nOfflatticeParticles() - 1) = y;
    m_particlePositions(2, nOfflatticeParticles() - 1) = z;

    m_F.resize(3, nOfflatticeParticles());

    m_localRates.resize(solver().length(), solver().width(), nOfflatticeParticles());

    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            m_localRates(x, y, nOfflatticeParticles() - 1) = calculateLocalRate(x, y, nOfflatticeParticles() - 1);
        }
    }

}

void OfflatticeMonteCarlo::initializeParticleMatrices(const uint nParticles, const double zMin)
{
    m_particlePositions.set_size(3, nParticles);
    m_F = zeros(3, nParticles);

    const double &h = solver().confiningSurfaceEvent().height();

    double x0;
    double y0;
    double z0;

    uint n = 0;
    while (n < nParticles)
    {
        do
        {
            double x0Raw = rng.uniform()*(solver().length() - 0.5);
            double y0Raw = rng.uniform()*(solver().width() - 0.5);
            z0 = zMin + rng.uniform()*(h - zMin);

            //it is rarely the case that the boundary transformations
            //touch a point inside the solver cube.
            x0 = solver().boundaryTransform(x0Raw, y0Raw, z0, 0);
            y0 = solver().boundaryTransform(x0, y0Raw, z0, 1);

        } while(solver().isBlockedPosition(x0, y0, z0));

        m_particlePositions(0, n) = x0;
        m_particlePositions(1, n) = y0;
        m_particlePositions(2, n) = z0;

        n++;
    }

    m_localRates.set_size(solver().length(), solver().width(), nOfflatticeParticles());

}

void OfflatticeMonteCarlo::scan(const uint n, const uint dim, const double dr, const uint maxSteps)
{
    const double &x0 = particlePositions(0, n);
    const double &y0 = particlePositions(1, n);
    const double &z0 = particlePositions(2, n);

    uint c = 0;

    while (solver().isBlockedPosition(x0, y0, z0) && c < maxSteps)
    {
        particlePositions(dim, n) = solver().boundaryTransform(particlePositions(0, n),
                                                               particlePositions(1, n),
                                                               particlePositions(2, n),
                                                               dr, dim);

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
            particlePositions(dim, n) = m_scanOriginalPositions(dim);
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

    lammpswriter writer(4, "cavitydiff", path);
    writer.setSystemSize(solver().length(), solver().width(), h, 0, 0, zMin);
    writer.initializeNewFile(frameNumber);

    if (nOfflatticeParticles() == 0)
    {
        writer << 0 << 0 << 0 << solver().height(0, 0);
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
               << z;
    }

    writer.finalize();
}

void OfflatticeMonteCarlo::clearDiffusingParticles()
{
    m_particlePositions.clear();
    m_F.clear();
    m_localRates.clear();
}

void OfflatticeMonteCarlo::calculateLocalRates()
{
    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            for (uint n = 0; n < nOfflatticeParticles(); ++n)
            {
                m_localRates(x, y, n) = calculateLocalRate(x, y, n);
            }
        }
    }
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

void OfflatticeMonteCarlo::setupInitialConditions()
{
    const double V = solver().volume();

    uint N = V*solver().concentration();

    const double zMin = solver().heights().min();
    initializeParticleMatrices(N, zMin);

    calculateLocalRates();
}



void OfflatticeMonteCarlo::execute()
{
}
