#include "offlatticemontecarlo.h"
#include "../confiningsurface/confiningsurface.h"
#include "../../solidonsolidreaction.h"

OfflatticeMonteCarlo::OfflatticeMonteCarlo(SolidOnSolidSolver &solver,
                                           const double D,
                                           const double dt) :
    Diffusion(solver, "OfflatticeMC"),
    m_D0(D),
    m_D(D*exp(-solver.gamma())),
    m_dt0(dt),
    m_dt(dt*exp(solver.gamma())/solver.area())
{

}

OfflatticeMonteCarlo::~OfflatticeMonteCarlo()
{

}

void OfflatticeMonteCarlo::dump(const uint frameNumber) const
{
    const double &h = solver().confiningSurfaceEvent().height();
    const int zMin = solver().heights().min();

    lammpswriter writer(4, "cavitydiff", "/tmp");
    writer.setSystemSize(solver().length(), solver().width(), h, 0, 0, zMin);
    writer.initializeNewFile(frameNumber);

    lammpswriter surfacewriter(5, "surfaces", "/tmp");
    surfacewriter.setSystemSize(solver().length(), solver().width(), h, 0, 0, zMin);
    surfacewriter.initializeNewFile(frameNumber);


    for (uint n = 0; n < nParticles(); ++n)
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

    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            surfacewriter << 1
                          << x
                          << y
                          << h
                          << 5;

            for (int zSurface = zMin; zSurface <= solver().height(x, y) - 1; ++zSurface)
            {
                surfacewriter << 2
                              << x
                              << y
                              << zSurface
                              << 6;
            }

            surfacewriter << 2
                          << x
                          << y
                          << solver().height(x, y)
                          << solver().nNeighbors(x, y);
        }
    }

    writer.finalize();
    surfacewriter.finalize();
}

void OfflatticeMonteCarlo::diffuse(const double dt)
{
    if (dt < 0 || fabs(dt) < 1E-15)
    {
        return;
    }

    double x1, y1, z1;
    double x0, y0, z0;

    for (uint n = 0; n < nParticles(); ++n)
    {
        x0 = m_particlePositions(0, n);
        y0 = m_particlePositions(1, n);
        z0 = m_particlePositions(2, n);

        BADAssBool(!solver().isBlockedPosition(x0, y0, z0));

        m_F(2, n) = solver().confiningSurfaceEvent().diffusionDrift(x0, y0, z0);

        do
        {
            x1 = x0 + sqrt(2*m_D*dt)*rng.normal() + m_D*m_F(0, n)*dt;
            y1 = y0 + sqrt(2*m_D*dt)*rng.normal() + m_D*m_F(1, n)*dt;
            z1 = z0 + sqrt(2*m_D*dt)*rng.normal() + m_D*m_F(2, n)*dt;

            if (solver().dim() == 1)
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

void OfflatticeMonteCarlo::removeParticle(const uint n)
{
    m_particlePositions.shed_col(n);
    m_F.shed_col(n);
    m_localProbabilities.shed_slice(n);
}

void OfflatticeMonteCarlo::insertParticle(const double x, const double y, const double z)
{
    m_particlePositions.resize(3, nParticles() + 1);

    m_particlePositions(0, nParticles() - 1) = x;
    m_particlePositions(1, nParticles() - 1) = y;
    m_particlePositions(2, nParticles() - 1) = z;

    m_F.resize(3, nParticles());
    m_localProbabilities.resize(solver().length(), solver().width(), nParticles());
}

double OfflatticeMonteCarlo::calculateTimeStep(const double initialCondition, bool calculateDissolutionRate)
{
    double timeStep = initialCondition;

    double totalDissolutionRate = 0;

    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            if (calculateDissolutionRate)
            {
                totalDissolutionRate += solver().reaction(x, y).calculateDissolutionRate();
            }

            else
            {
                totalDissolutionRate += solver().reaction(x, y).dissolutionRate();
            }
        }
    }

    double eps;
    uint i = 0;

    cout << "start: " << timeStep << endl;

    do
    {
        double totalDepositionRate = 0;
        for (uint x = 0; x < solver().length(); ++x)
        {
            for (uint y = 0; y < solver().width(); ++y)
            {
                totalDepositionRate += depositionRate(x, y, timeStep);
            }
        }

        double newTimeStep = solver().nextRandomLogNumber()/(totalDissolutionRate + totalDepositionRate);

        cout << i << " " << newTimeStep << endl;

        i++;

        eps = fabs(timeStep - newTimeStep);

        timeStep = newTimeStep;

    } while (eps > m_eps);

    double totalDepositionRate = 0;
    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            totalDepositionRate += depositionRate(x, y, timeStep);
        }
    }

    double newTimeStep = solver().nextRandomLogNumber()/(totalDissolutionRate + totalDepositionRate);

    cout << fabs(newTimeStep - timeStep) << endl;
    cout << "---" << endl;

    return timeStep;
}

double OfflatticeMonteCarlo::depositionRate(const uint x, const uint y, double timeStep) const
{
    double P = 0;

    const double sigmaSquared = 2*solver().dim()*m_D*timeStep;

    for (uint n = 0; n < nParticles(); ++n)
    {

        const double Pn = calculateLocalProbability(x, y, n, timeStep);

        P += Pn;
    }

    double geometric = (6 - solver().nNeighbors(x, y))/(4*datum::pi);
    return P*exp(2*solver().alpha() - solver().gamma())/sqrt(2*datum::pi*sigmaSquared)*geometric;
}

double OfflatticeMonteCarlo::calculateLocalProbability(const uint x,
                                                       const uint y,
                                                       const uint n,
                                                       const double timeStep) const
{

    const int z = solver().height(x, y) + 1;

    const double sigmaSquared = 2*solver().dim()*m_D*timeStep;

    const double dxSquared = pow((double)x - m_particlePositions(0, n), 2);
    const double dySquared = pow((double)y - m_particlePositions(1, n), 2);
    const double dzSquared = pow((double)z - m_particlePositions(2, n), 2);

    const double dr2 = dxSquared + dySquared + dzSquared;

    BADAss(dr2, !=, 0, "lols", [&] () {
        cout << m_particlePositions << endl;
    });

    return exp(-dr2/(2*sigmaSquared))/dr2;
}



void OfflatticeMonteCarlo::initialize()
{
    m_accepted = 0;
    m_trials = 0;
}

void OfflatticeMonteCarlo::execute()
{
    m_currentTimeStep = calculateTimeStep(m_currentTimeStep);

    if (cycle() % 100 == 0)
    {
        dump(cycle()/100);
    }
}

void OfflatticeMonteCarlo::reset()
{
    //Has not yet been updated: Represents time spent in current state.
    const double &dtFull = solver().currentTimeStep();

    //    BADAssClose(dtFull, m_currentTimeStep, 1E-3); //DERP

    const uint N = dtFull/m_dt;

    for (uint i = 0; i < N; ++i)
    {
        diffuse(m_dt);
    }

    diffuse(dtFull - N*m_dt);
}


void OfflatticeMonteCarlo::setupInitialConditions()
{
    const double V = solver().volume();
    const double &h = solver().confiningSurfaceEvent().height();

    uint N = V*exp(solver().gamma()-2*solver().alpha());
    N = 10000;

    m_particlePositions.set_size(3, N);

    double zMin = solver().heights().min();

    double x0;
    double y0;
    double z0;

    uint n = 0;
    while (n < N)
    {
        do
        {
            x0 = rng.uniform()*(solver().length() - 1);
            y0 = rng.uniform()*(solver().width() - 1);
            z0 = zMin + rng.uniform()*(h - zMin);

        } while(solver().isBlockedPosition(x0, y0, z0));

        m_particlePositions(0, n) = x0;
        m_particlePositions(1, n) = y0;
        m_particlePositions(2, n) = z0;

        n++;
    }

    m_F = zeros(3, N);

    m_localProbabilities.set_size(solver().length(), solver().width(), nParticles());

    m_currentTimeStep = calculateTimeStep(1.0, true);

}

void OfflatticeMonteCarlo::registerHeightChange(const uint x, const uint y, const int delta)
{
    if (!hasStarted())
    {
        return;
    }

    //Remove particle based on the probability of it being the deposited
    if (delta == 1)
    {
        vector<double> localRatesForSite(nParticles());

        double Rtot = 0;
        for (uint n = 0; n < nParticles(); ++n)
        {
            m_localProbabilities(x, y, n) = calculateLocalProbability(x, y, n);
            BADAssClose(m_localProbabilities(x, y, n), calculateLocalProbability(x, y, n), 1E-5);

            localRatesForSite.at(n) = m_localProbabilities(x, y, n);
            Rtot += localRatesForSite.at(n);
        }

        double R = rng.uniform()*Rtot;

        uint N = binarySearchForInterval(R, localRatesForSite);

        removeParticle(N);
    }

    //Add a particle on a unit sphere around the site;
    else
    {
        double x1, y1, z1;

        const int &z = solver().height(x, y) + 1;

        do
        {
            double theta = rng.uniform()*datum::pi;
            double phi = rng.uniform()*datum::pi*2;

            x1 = (double)x + sin(theta)*cos(phi);
            y1 = (double)y + sin(theta)*sin(phi);
            z1 = (double)z + cos(theta);

        } while (solver().isBlockedPosition(x1, y1, z1));

        insertParticle(x1, y1, z1);
    }

    const double &h = solver().confiningSurfaceEvent().height();
    double zMin = solver().heights().min();

    for (uint n = 0; n < nParticles(); ++n)
    {
        double &x0 = m_particlePositions(0, n);
        double &y0 = m_particlePositions(1, n);
        double &z0 = m_particlePositions(2, n);

        while (solver().isBlockedPosition(x0, y0, z0))
        {
            x0 = rng.uniform()*solver().length();

            if (solver().dim() != 1)
            {
                y0 = rng.uniform()*solver().width();
            }

            z0 = zMin + rng.uniform()*(h - zMin);
        }
    }

    m_currentTimeStep = calculateTimeStep(m_currentTimeStep);

    DiffusionDeposition *r;
    for (uint _x = 0; _x < solver().length(); ++_x)
    {
        for (uint _y = 0; _y < solver().width(); ++_y)
        {
            r = &solver().reaction(_x, _y);
            r->setDepositionRate(r->calculateDepositionRate());
        }
    }
}

