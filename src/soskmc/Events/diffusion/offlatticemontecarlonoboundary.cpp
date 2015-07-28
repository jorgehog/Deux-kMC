#include "offlatticemontecarlonoboundary.h"
#include "../confiningsurface/confiningsurface.h"
#include "../../dissolutiondeposition.h"

OfflatticeMonteCarloNoBoundary::OfflatticeMonteCarloNoBoundary(SOSSolver &solver,
                                                               const double dt) :
    OfflatticeMonteCarlo(solver, dt, "OfflatticeMC")
{

}

OfflatticeMonteCarloNoBoundary::~OfflatticeMonteCarloNoBoundary()
{

}


double OfflatticeMonteCarloNoBoundary::calculateTimeStep(const double initialCondition, bool calculateDissolutionRate)
{
    double timeStep = initialCondition;

    double totalDissolutionRate = 0;

    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            if (calculateDissolutionRate)
            {
                totalDissolutionRate += solver().surfaceReaction(x, y).calculateDissolutionRate();
            }

            else
            {
                totalDissolutionRate += solver().surfaceReaction(x, y).dissolutionRate();
            }
        }
    }

    double eps;
    uint i = 0;

    vector<double> prevTimeSteps;
    double totalDepositionRate;
    double _depositionRate;
    do
    {
        totalDepositionRate = 0;
        for (uint x = 0; x < solver().length(); ++x)
        {
            for (uint y = 0; y < solver().width(); ++y)
            {
                _depositionRate = depositionRate(x, y, timeStep);
                totalDepositionRate += _depositionRate;
            }
        }

        double newTimeStep = solver().nextRandomLogNumber()/(totalDissolutionRate + totalDepositionRate);

        cout << "\r" << std::setw(3) << i << " " << std::setprecision(6) << std::fixed << newTimeStep;
        cout.flush();

        i++;

        eps = fabs(timeStep - newTimeStep);

        timeStep = newTimeStep;

        //checks the new time step against previous time steps
        //if it is identical, then the iteration failed.
        for (const double &prevTimeStep : prevTimeSteps)
        {
            if (fabs(prevTimeStep - timeStep) < 1E-15)
            {
                cout << "Fixed point iteration failed." << prevTimeStep << " " << timeStep << endl;

                for (const double &prevTimeStep2 : prevTimeSteps)
                {
                    cout << prevTimeStep2 << " ";
                }
                cout << endl;
                sleep(2);
                return solver().nextRandomLogNumber()/totalDissolutionRate;
            }
        }


        //cycle previous time steps back one index and add the new one to the end.
        if (!prevTimeSteps.empty())
        {
            for (uint j = 0; j < prevTimeSteps.size() - 1; ++j)
            {
                prevTimeSteps.at(j) = prevTimeSteps.at(j + 1);
            }
            prevTimeSteps.back() = timeStep;
        }

        //every 100 iteration cycles we check against one extra time step
        //to be able to locate larger loops, i.e. a-b-c-d-e-f-d-e-f-d-e-f-.. etc.
        //and not only let's say a-b-c-b-c-b-c..
        if (i % 20 == 0)
        {
            prevTimeSteps.push_back(timeStep);
        }

    } while (eps > m_eps);

    //    double totalDepositionRate = 0;
    //    for (uint x = 0; x < solver().length(); ++x)
    //    {
    //        for (uint y = 0; y < solver().width(); ++y)
    //        {
    //            totalDepositionRate += depositionRate(x, y, timeStep);
    //        }
    //    }

    //    double newTimeStep = solver().nextRandomLogNumber()/(totalDissolutionRate + totalDepositionRate);

    cout << "  " << eps << " " << timeStep/initialCondition << " --- " << totalDepositionRate << " " << totalDissolutionRate << " " << totalDepositionRate/totalDissolutionRate << endl;

    return timeStep;
}

double OfflatticeMonteCarloNoBoundary::depositionRate(const uint x, const uint y, double timeStep) const
{
    double P = 0;

    for (uint n = 0; n < nOfflatticeParticles(); ++n)
    {

        const double Pn = calculateLocalProbability(x, y, n, timeStep);

        P += Pn;
    }

    return P/solver().concentration();
}

double OfflatticeMonteCarloNoBoundary::calculateLocalProbability(const uint x,
                                                                 const uint y,
                                                                 const uint n,
                                                                 const double timeStep) const
{

    const int z = solver().height(x, y) + 1;

    const double sigmaSquared = 4*D()*timeStep;

    const double dxSquared = pow((double)x - particlePositions(0, n), 2);
    const double dySquared = pow((double)y - particlePositions(1, n), 2);
    const double dzSquared = pow((double)z - particlePositions(2, n), 2);

    const double dr2 = dxSquared + dySquared + dzSquared;
    //    const double r = sqrt(dr2);

    BADAss(dr2, !=, 0, "lols", [&] () {
        cout << particlePositions() << endl;
    });

    if (solver().dim() == 2)
    {
        BADAssBreak("Not supported in 2D.");
        cout << "2D sucks" << endl;
        exit(1);
    }

    //This takes into account that the particle can deposit at any time between t and t + dt, not only at t + dt.
    //    const double N = 4*m_D*r;
    //    const double N = r;
    //    const double N = 1;
    //    return std::erfc(r/sqrt(sigmaSquared))/N;

    return exp(-dr2/(2*sigmaSquared))/sqrt(2*datum::pi*sigmaSquared);
}


void OfflatticeMonteCarloNoBoundary::execute()
{
    m_currentTimeStep = calculateTimeStep(m_currentTimeStep);

    const uint THRESH = 1;
    if (cycle() % THRESH == 0)
    {
        dump(cycle()/THRESH);
    }
}

void OfflatticeMonteCarloNoBoundary::reset()
{
    //Has not yet been updated: Represents time spent in current state.
    const double &dtFull = solver().currentTimeStep();

    //    BADAssClose(dtFull, m_currentTimeStep, 1E-3); //DERP

    const uint N = dtFull/dt();

    for (uint i = 0; i < N; ++i)
    {
        diffuse(dt());
    }

    diffuse(dtFull - N*dt());
}


void OfflatticeMonteCarloNoBoundary::setupInitialConditions()
{
    const double V = solver().volume();

    uint N = V*solver().concentration();

    const double zMin = solver().heights().min();
    initializeParticleMatrices(N, zMin);

    m_localProbabilities.set_size(solver().length(), solver().width(), nOfflatticeParticles());

    m_currentTimeStep = calculateTimeStep(1.0, true);

}

void OfflatticeMonteCarloNoBoundary::registerHeightChange(const uint x, const uint y, const int delta)
{
    if (!hasStarted())
    {
        return;
    }

    //Remove particle based on the probability of it being the deposited
    if (delta == 1)
    {
        vector<double> localRatesForSite(nOfflatticeParticles());

        double Rtot = 0;
        for (uint n = 0; n < nOfflatticeParticles(); ++n)
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

    for (uint n = 0; n < nOfflatticeParticles(); ++n)
    {
        double &x0 = particlePositions(0, n);
        double &y0 = particlePositions(1, n);
        double &z0 = particlePositions(2, n);

        while (solver().isBlockedPosition(x0, y0, z0))
        {
            x0 = rng.uniform()*solver().length();

            if (solver().surfaceDim() != 1)
            {
                y0 = rng.uniform()*solver().width();
            }

            z0 = zMin + rng.uniform()*(h - zMin);
        }
    }

    m_currentTimeStep = calculateTimeStep(m_currentTimeStep);

    DissolutionDeposition *r;
    for (uint _x = 0; _x < solver().length(); ++_x)
    {
        for (uint _y = 0; _y < solver().width(); ++_y)
        {
            r = &solver().surfaceReaction(_x, _y);
            r->setDepositionRate(r->calculateDepositionRate());
        }
    }
}

void OfflatticeMonteCarloNoBoundary::executeDiffusionReaction(SOSDiffusionReaction *reaction,
                                                              const uint x, const uint y, const int z)
{
    (void) reaction;
    (void) x;
    (void) y;
    (void) z;

    BADAssBreak("Invalid diffusion scheme for diffusion reactions");
}

bool OfflatticeMonteCarloNoBoundary::isBlockedPosition(const uint x, const uint y, const int z) const
{
    (void) x;
    (void) y;
    (void) z;

    return false;
}



void OfflatticeMonteCarloNoBoundary::onInsertParticle(const double x, const double y, const double z)
{
    (void) x;
    (void) y;
    (void) z;

    m_localProbabilities.resize(solver().length(), solver().width(), nOfflatticeParticles());
}

void OfflatticeMonteCarloNoBoundary::onRemoveParticle(const uint n)
{
    m_localProbabilities.shed_slice(n);
}

