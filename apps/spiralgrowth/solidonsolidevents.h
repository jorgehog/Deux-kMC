#pragma once

#include <kMC>
#include <utils.h>
#include "solidonsolidsolver.h"

using namespace ignis;

class SolidOnSolidEvent : public ignis::LatticeEvent
{
public:
    using LatticeEvent::LatticeEvent;
    ~SolidOnSolidEvent();

    const SolidOnSolidSolver *solver() const
    {
        return dependency<SolidOnSolidSolver>("KMCSolver");
    }

};

class SurfaceSize : public SolidOnSolidEvent
{
public:

    SurfaceSize() :
        SolidOnSolidEvent("SurfaceSize", "l0", true, true)
    {

    }

    void initialize()
    {
        m_sum = 0;
        m_T0 = solver()->currentTime() - solver()->currentTimeStep();
    }

    void execute();

    const double &localValue() const
    {
        return m_localValue;
    }

private:

    double m_sum;
    double m_localValue;

    double m_T0;

};


class AverageHeight : public SolidOnSolidEvent
{
public:

    AverageHeight()
        : SolidOnSolidEvent("AverageHeight", "l0", true, true)
    {

    }

    void execute();
};


class DumpHeights3D : public SolidOnSolidEvent
{
public:

    DumpHeights3D(const string path = "/tmp") :
        SolidOnSolidEvent("DumpHeights3D"),
        m_writer(4, "SOS", path)
    {

    }

    void initialize();

    void execute();

private:

    const string m_filename;

    lammpswriter m_writer;

    uint m_L;
    uint m_W;
    uint m_N;

};


class NNeighbors : public SolidOnSolidEvent
{
public:

    NNeighbors() :
        SolidOnSolidEvent("nNeighbors", "", true, true)
    {

    }

    const double &localValue() const
    {
        return m_localValue;
    }


    // Event interface
public:
    void execute();

private:

    double m_localValue;

};


class EqMu : public SolidOnSolidEvent
{
public:

    EqMu() :
        SolidOnSolidEvent("EqMu", "", true, true),
        m_accuNeighbours(0),
        m_accuDissolutionRate(0),
        m_totalTime(0)
    {

    }

    ~EqMu()
    {

    }

    void initialize();

    void execute();

    void reset();

    void restart();

    double dMu() const
    {
        return log(m_dMu);
    }

private:

    double m_dMu;

    double m_accuNeighbours;
    double m_accuDissolutionRate;

    double m_totalTime;

    void update();

    void resetCounters()
    {
        m_totalTime = 0;

        m_accuNeighbours = 0;
        m_accuDissolutionRate = 0;
    }


};

class Equilibriater : public SolidOnSolidEvent
{
public:

    Equilibriater(SolidOnSolidSolver &mutexSolver,
                  EqMu &eqMuEvent,
                  const uint nRounds = 100,
                  const uint N = 1000) :
        SolidOnSolidEvent("ConcEquilibriator", "", true, true),
        m_mutexSolver(mutexSolver),
        m_eqMuEvent(eqMuEvent),
        m_nRounds(nRounds),
        m_N(N),
        m_doAverage(false),
        m_averageMu(0),
        m_averageMu2Sum(0),
        m_averageMuCount(0),
        m_finalized(false)
    {
        setDependency(m_eqMuEvent);
    }

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

    SolidOnSolidSolver &m_mutexSolver;

    EqMu &m_eqMuEvent;

    const uint m_nRounds;
    const uint m_N;

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

