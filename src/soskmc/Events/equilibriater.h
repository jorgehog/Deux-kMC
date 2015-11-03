#pragma once

#include "sosevent.h"

class EqGamma;

class Equilibriater : public SOSEvent
{
public:

    Equilibriater(SOSSolver &mutexSolver,
                  EqGamma &eqMuEvent,
                  const uint nSamplesGammaEq = 100,
                  const uint nSamplesGamma = 1000);

    void initialize();

    void execute();

    const double &averageGamma() const
    {
        BADAssBool(m_finalized);

        return m_averageGamma;
    }

    const double &error() const
    {
        BADAssBool(m_finalized);

        return m_error;
    }

    const bool &finalized() const
    {
        return m_finalized;
    }

    void finalizeAverages();

private:

    SOSSolver &m_mutexSolver;

    EqGamma &m_eqGammaEvent;

    const uint m_nSamplesGammaEq;
    const uint m_nSamplesGamma;

    bool m_doAverage;
    double m_averageGamma;
    double m_averageGamma2Sum;
    double m_error;
    uint m_averageGammaCount;

    bool m_finalized;

    uint m_counter;

    void initiateNextConcentrationLevel(const double shift);

    vector<double> m_shifts;
    vector<double> m_values;

    double m_prevShift;

};

