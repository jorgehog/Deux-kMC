#include "equilibriater.h"

#include "eqmu.h"


Equilibriater::Equilibriater(SolidOnSolidSolver &mutexSolver,
                             EqMu &eqMuEvent,
                             const uint nSamplesMuEq,
                             const uint nSamplesMu) :
    SolidOnSolidEvent(mutexSolver,
                      "ConcEquilibriator",
                      "",
                      true,
                      true),
    m_mutexSolver(mutexSolver),
    m_eqMuEvent(eqMuEvent),
    m_nSamplesMuEq(nSamplesMuEq),
    m_nSamplesMu(nSamplesMu),
    m_doAverage(false),
    m_averageMu(0),
    m_averageMu2Sum(0),
    m_averageMuCount(0),
    m_finalized(false)
{
    setDependency(m_eqMuEvent);
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
    if (m_counter >= m_nSamplesMu)
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


    if (m_averageMuCount == m_nSamplesMuEq)
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

    m_finalized = true;

}
