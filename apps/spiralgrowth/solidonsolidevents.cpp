#include "solidonsolidevents.h"

#include "solidonsolidreaction.h"

SolidOnSolidEvent::~SolidOnSolidEvent()
{

}

void SurfaceSize::execute()
{
    m_localValue = 0;

    for (uint x = 0; x < solver()->length(); ++x)
    {
        for (uint y = 0; y < solver()->width(); ++y)
        {
            m_localValue += 0.5*abs(solver()->height(x, y) - solver()->height(solver()->rightSite(x), y));
            m_localValue += 0.5*abs(solver()->height(x, y) - solver()->height(x, solver()->topSite(y)));
        }
    }

    m_localValue /= (solver()->length()*solver()->width());

    m_sum += solver()->currentTimeStep()*m_localValue;

    setValue(m_sum/(solver()->currentTime() - m_T0));
}


void DumpHeights3D::initialize()
{
    m_L = solver()->length();
    m_W = solver()->width();
    m_N = m_L*m_W;
}

void DumpHeights3D::execute()
{
    if (cycle() % 100 != 0)
    {
        return;
    }

    const imat &heights = solver()->heights();

    const int max = heights.max();
    const int min = heights.min();
    const uint maxSpan = max - min;

    m_writer.setSystemSize(m_L, m_W, maxSpan, 0, 0, 0);

    m_writer.initializeNewFile(cycle()/100);

    for (uint x = 0; x < m_L; ++x)
    {
        for (uint y = 0; y < m_W; ++y)
        {
            uint zSpan = heights(x, y) - min;

            for (uint z = 0; z < zSpan; ++z)
            {
                m_writer << x
                         << y
                         << z
                         << 0;
            }

            m_writer << x
                     << y
                     << zSpan
                     << solver()->nNeighbors(x, y);
        }
    }

    m_writer.finalize();

}


void AverageHeight::execute()
{
    setValue(accu(solver()->heights())/(double)solver()->area());
}



void EqMu::initialize()
{
    resetCounters();

    update();
}

void EqMu::execute()
{
    setValue(dMu());
}

void EqMu::reset()
{
    update();
}

void EqMu::restart()
{
    resetCounters();
}

void EqMu::update()
{
    double localDissolutionRate = 0;

    for (uint x = 0; x < solver()->length(); ++x)
    {
        for (uint y = 0; y < solver()->width(); ++y)
        {
            localDissolutionRate += solver()->reaction(x, y).diffusionRate();
        }
    }

    localDissolutionRate /= solver()->area();

    const double &dt = solver()->currentTimeStep();

    m_totalTime += dt;

    m_accuDissolutionRate += dt*localDissolutionRate;

    double avgDissolutionRate = m_accuDissolutionRate/m_totalTime;

    double inverseKStar;

    if (solver()->shadowing())
    {
        m_accuNeighbours += dt*dependency<NNeighbors>("nNeighbors")->localValue();
        const double &avgNeighbors = m_accuNeighbours/m_totalTime;

        inverseKStar = 1./solver()->shadowScale(avgNeighbors);
    }
    else
    {
        inverseKStar = 1.;
    }

    m_dMu = avgDissolutionRate*inverseKStar;
}


void Equilibriater::initialize()
{
    m_counter = 0;

    m_prevShift = 0;

    m_doAverage = false;

    m_averageMu = 0;

    m_averageMu2Sum = 0;

    m_averageMuCount = 0;

    m_finalized = false;

}

void Equilibriater::execute()
{
    double shift = m_eqMuEvent.value();

    m_counter++;
    if (m_counter >= m_N)
    {
        initiateNextConcentrationLevel(shift);
        m_counter = 0;
    }

    setValue(m_mutexSolver.mu());

}


void Equilibriater::initiateNextConcentrationLevel(const double shift)
{

    double newMu = m_mutexSolver.mu() + shift;

    //Cannot use 'else' because this could occur even when the prev test goes through
    if (m_doAverage)
    {
        m_averageMu += newMu;
        m_averageMu2Sum += newMu*newMu;
        m_averageMuCount++;
    }

    m_shifts.push_back(shift);
    m_values.push_back(m_mutexSolver.mu());

    m_mutexSolver.setMu(newMu);

    conv_to<vec>::from(m_shifts).eval().save("/tmp/shifts.arma");
    conv_to<vec>::from(m_values).eval().save("/tmp/values.arma");


    if (m_averageMuCount == m_nRounds)
    {
        terminateLoop("Concentration converged");

        finalizeAverages();
    }

    if (!m_doAverage)
    {
        if (m_prevShift != 0 && shift != 0)
        {
            if (m_prevShift/shift < 0)
            {
                m_doAverage = true;
            }
        }
    }

    m_prevShift = shift;

    m_eqMuEvent.restart();

    for (uint n = 0; n < m_mutexSolver.area(); ++n)
    {
        m_mutexSolver.getReaction(n)->calculateRate();
    }

}

void Equilibriater::finalizeAverages()
{
    if (m_finalized)
    {
        return;
    }

    for (uint i = 0; i < m_shifts.size(); ++i)
    {
        cout << m_shifts.at(i) << " " << m_values.at(i) << endl;
    }

    const uint &N = m_averageMuCount;

    cout << "N = " << N << endl;

    m_averageMu /= N;

    m_error = sqrt(1.0/(N - 1)*(m_averageMu2Sum - m_averageMu*m_averageMu*N));

    cout << "Average = " << m_averageMu << " error = " << m_error << endl;

    m_finalized = true;

    //Update rates?
}



void NNeighbors::execute()
{
}
