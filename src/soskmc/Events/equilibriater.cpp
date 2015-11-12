#include "sossolver.h"

#include "equilibriater.h"

#include "eqmu.h"



Equilibriater::Equilibriater(SOSSolver &mutexSolver,
                             EqGamma &eqMuEvent,
                             const uint nSamplesGammaEq,
                             const uint nSamplesGamma) :
    SOSEvent(mutexSolver,
                      "ConcEquilibriator",
                      "",
                      true,
                      true),
    m_mutexSolver(mutexSolver),
    m_eqGammaEvent(eqMuEvent),
    m_nSamplesGammaEq(nSamplesGammaEq),
    m_nSamplesGamma(nSamplesGamma),
    m_doAverage(false),
    m_averageGamma(0),
    m_averageGamma2Sum(0),
    m_averageGammaCount(0),
    m_finalized(false)
{
    setDependency(m_eqGammaEvent);
}

void Equilibriater::initialize()
{
    m_counter = 0;

    m_prevShift = 0;

    m_doAverage = false;

    m_averageGamma = 0;

    m_averageGamma2Sum = 0;

    m_averageGammaCount = 0;

    m_finalized = false;

}

void Equilibriater::execute()
{
    double shift = m_eqGammaEvent.value();

    m_counter++;
    if (m_counter >= m_nSamplesGamma)
    {
        initiateNextConcentrationLevel(shift);
        m_counter = 0;
    }

    setValue(m_mutexSolver.gamma());

}


void Equilibriater::initiateNextConcentrationLevel(const double shift)
{
    double newGamma = m_mutexSolver.gamma() + shift;

    //Cannot use 'else' because this could occur even when the prev test goes through
    if (m_doAverage)
    {
        m_averageGamma += newGamma;
        m_averageGamma2Sum += newGamma*newGamma;
        m_averageGammaCount++;
    }

    m_shifts.push_back(shift);
    m_values.push_back(m_mutexSolver.gamma());

    m_mutexSolver.setGamma(newGamma);

    conv_to<vec>::from(m_shifts).eval().save("/tmp/shifts.arma");
    conv_to<vec>::from(m_values).eval().save("/tmp/values.arma");


    if (m_averageGammaCount == m_nSamplesGammaEq)
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

    m_eqGammaEvent.restart();

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

    const uint &N = m_averageGammaCount;

    m_averageGamma /= N;

    m_error = sqrt(1.0/(N - 1)*(m_averageGamma2Sum - m_averageGamma*m_averageGamma*N));

    m_finalized = true;

}
