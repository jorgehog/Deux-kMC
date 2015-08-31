#pragma once

#include "sosevent.h"

class EqMu;

class Equilibriater : public SOSEvent
{
public:

    Equilibriater(SOSSolver &mutexSolver,
                  EqMu &eqMuEvent,
                  const uint nSamplesMuEq = 100,
                  const uint nSamplesMu = 1000);

    void initialize();

    void execute();

    const double &averageMu() const
    {
        BADAssBool(m_finalized);

        return m_averageMu;
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

    EqMu &m_eqMuEvent;

    const uint m_nSamplesMuEq;
    const uint m_nSamplesMu;

    bool m_doAverage;
    double m_averageMu;
    double m_averageMu2Sum;
    double m_error;
    uint m_averageMuCount;

    bool m_finalized;

    uint m_counter;

    void initiateNextConcentrationLevel(const double shift);

    vector<double> m_shifts;
    vector<double> m_values;

    double m_prevShift;

};

