#include "fixedpointtimestepping.h"

#include "Events/confiningsurface/confiningsurface.h"
#include "dissolutiondeposition.h"

#include "sossolver.h"

FixedPointTimeStepping::FixedPointTimeStepping(SOSSolver &solver, const double maxdt) :
    Diffusion(solver, "FixedPoint"),
    OfflatticeMonteCarlo(solver, maxdt, 0),
    m_mutexSolver(solver)
{

}

void FixedPointTimeStepping::calculateTimeStep(bool calculateDissolutionRate)
{
    double totalDissolutionRate = 0;

    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            if (calculateDissolutionRate)
            {
                totalDissolutionRate += solver().surfaceReaction(x, y).calculateEscapeRate();
            }

            else
            {
                totalDissolutionRate += solver().surfaceReaction(x, y).escapeRate();
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
                _depositionRate = depositionRate(x, y);
                totalDepositionRate += _depositionRate;
            }
        }

        double newTimeStep = solver().nextRandomLogNumber()/(totalDissolutionRate + totalDepositionRate);

        cout << "\r" << std::setw(3) << i << " " << std::setprecision(6) << std::fixed << newTimeStep;
        cout.flush();

        i++;

        eps = fabs(m_currentTimeStep - newTimeStep);

        m_currentTimeStep = newTimeStep;

        //checks the new time step against previous time steps
        //if it is identical, then the iteration failed.
        for (const double &prevTimeStep : prevTimeSteps)
        {
            if (fabs(prevTimeStep - m_currentTimeStep) < 1E-15)
            {
                cout << "Fixed point iteration failed." << prevTimeStep << " " << m_currentTimeStep << endl;

                for (const double &prevTimeStep2 : prevTimeSteps)
                {
                    cout << prevTimeStep2 << " ";
                }
                cout << endl;
                sleep(2);
                return;
            }
        }


        //cycle previous time steps back one index and add the new one to the end.
        if (!prevTimeSteps.empty())
        {
            for (uint j = 0; j < prevTimeSteps.size() - 1; ++j)
            {
                prevTimeSteps.at(j) = prevTimeSteps.at(j + 1);
            }
            prevTimeSteps.back() = m_currentTimeStep;
        }

        //every 100 iteration cycles we check against one extra time step
        //to be able to locate larger loops, i.e. a-b-c-d-e-f-d-e-f-d-e-f-.. etc.
        //and not only let's say a-b-c-b-c-b-c..
        if (i % 20 == 0)
        {
            prevTimeSteps.push_back(m_currentTimeStep);
        }

        calculateLocalRatesAndUpdateDepositionRates();

    } while (eps > m_eps);


    cout << "  " << eps << " --- " << totalDepositionRate << " " << totalDissolutionRate << " " << totalDepositionRate/totalDissolutionRate << endl;
}

double FixedPointTimeStepping::calculateLocalRateOverD(const uint x,
                                                  const uint y,
                                                  const uint n,
                                                  const double timeStep) const
{

    const int z = solver().height(x, y) + 1;

    const double sigmaSquared = 4*DScaled()*timeStep;

    const double dxSquared = pow((double)x - particlePositions(0, n), 2);
    const double dySquared = pow((double)y - particlePositions(1, n), 2);
    const double dzSquared = pow((double)z - particlePositions(2, n), 2);

    const double dr2 = dxSquared + dySquared + dzSquared;
    const double r = sqrt(dr2);

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
    const double N = 4*DScaled()*r;
    //    const double N = r;
    //    const double N = 1;
    return std::erfc(r/sqrt(sigmaSquared))/N/timeStep;

    //    return exp(-dr2/(2*sigmaSquared))/sqrt(2*datum::pi*sigmaSquared);
}


void FixedPointTimeStepping::execute()
{
    calculateTimeStep();
}

void FixedPointTimeStepping::initializeObserver(const Subjects &subject)
{
    OfflatticeMonteCarlo::initializeObserver(subject);

    m_currentTimeStep = 1.0;
    calculateTimeStep(true);

}

void FixedPointTimeStepping::notifyObserver(const Subjects &subject)
{
    OfflatticeMonteCarlo::notifyObserver(subject);

    calculateTimeStep();

    SurfaceReaction *r;
    for (uint _x = 0; _x < solver().length(); ++_x)
    {
        for (uint _y = 0; _y < solver().width(); ++_y)
        {
            r = &m_mutexSolver.surfaceReaction(_x, _y);
            r->setDepositionRate(r->calculateDepositionRate());
        }
    }
}

double FixedPointTimeStepping::calculateLocalRateOverD(const uint x, const uint y, const uint n) const
{
    return calculateLocalRateOverD(x, y, n, m_currentTimeStep);
}



void FixedPointTimeStepping::calculateLocalRatesAndUpdateDepositionRates()
{
}
