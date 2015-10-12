#pragma once

#include "sosevent.h"

class Time : public SOSEvent
{
public:
    Time(const SOSSolver &solver) :
        SOSEvent(solver, "Time", "", false, true)
    {

    }

    void execute()
    {
        setValue(solver().currentTime());
    }
};

class GrowthSpeed : public SOSEvent
{
public:

    GrowthSpeed(const SOSSolver &solver) :
        SOSEvent(solver, "GrowthSpeed", "", true)
    {

    }

    void initialize();

    void execute();

private:

    double m_T0;
    double m_h0;
};

class SurfaceSize : public SOSEvent
{
public:

    SurfaceSize(const SOSSolver &solver) :
        SOSEvent(solver, "SurfaceSize", "l0", true, true),
        m_relativeHeightSums(solver.length(), solver.width())
    {

    }

    void initialize();

    double timeAverage() const
    {
        return m_sum/(solver().currentTime() - m_T0);
    }

    double relativeHeightSum(const uint x, const uint y) const;

    double getLocalValue() const
    {
        return m_localValue/solver().area();
    }

    void updateRelativeHeight(const uint x, const uint y);

    double bruteForceValue() const;

    void execute();

    void reset();

private:

    double m_sum;
    double m_localValue;

    mat m_relativeHeightSums;

    double m_T0;

};

class SurfaceVariance : public SOSEvent
{
public:

    SurfaceVariance(const SOSSolver &solver) :
        SOSEvent(solver, "SurfaceVariance", "", true, true)
    {

    }

    void initialize()
    {
        m_s2 = 0;
        m_T0 = solver().currentTime() - solver().currentTimeStep();
    }

    void execute();

private:

    double m_s2;

    double m_T0;
};

class HeightRMS : public SOSEvent
{
public:
    HeightRMS(const SOSSolver &solver) :
        SOSEvent(solver, "HeightRMS", "", true, true)
    {

    }

    void execute();
};


class AverageHeight : public SOSEvent
{
public:

    AverageHeight(const SOSSolver &solver) :
        SOSEvent(solver, "AverageHeight", "l0", true, true)
    {

    }

    void initialize();

    void execute();

    double getValue() const;
};


class DumpHeightSlice : public SOSEvent
{
public:
    DumpHeightSlice(const SOSSolver &solver,
                    const uint slicePosition=0,
                    const uint axis=0,
                    const string path = "/tmp",
                    const uint nCyclesperOutput = 1000) :
        SOSEvent(solver, "DumpHeightSlice"),
        m_slicePosition(slicePosition),
        m_axis(axis),
        m_path(path),
        m_nCyclesPerOutput(nCyclesperOutput)
    {

    }

    void initialize();

    void execute();

private:

    const uint m_slicePosition;
    const uint m_axis;

    const string m_path;
    const uint m_nCyclesPerOutput;

    ivec m_heights;

};

class DumpHeights3D : public SOSEvent
{
public:

    DumpHeights3D(const SOSSolver &solver,
                  const string path = "/tmp") :
        SOSEvent(solver, "DumpHeights3D"),
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


class NNeighbors : public SOSEvent
{
public:

    NNeighbors(const SOSSolver &solver) :
        SOSEvent(solver, "nNeighbors", "", true, true),
        m_nNeighbors(solver.length(), solver.width())
    {

    }

    void updateNNeighbors(const uint x, const uint y);

    double bruteForceValue() const;

    // Event interface
public:
    void initialize();
    void execute();
    void reset();

private:

    double m_localValue;

    //need an additional neighbor matrix to store the old values.
    umat m_nNeighbors;

};


class RateChecker : public LatticeEvent
{

public:

    RateChecker(const KMCSolver &solver);

    void execute() {}

    void reset();

private:

    const KMCSolver &m_solver;

};


class AutoCorrelationHeight : public SOSEvent
{

public:

    AutoCorrelationHeight(const SOSSolver &solver,
                          const uint xSpan,
                          const uint ySpan);

    mat autoCorrelation() const;


private:

    const uint m_xSpan;
    const uint m_ySpan;

    mat m_autoCorrelationQuadrant;
    mat m_autoCorrelationSubQuadrant;

    // Event interface
public:
    void execute();

    void initialize();
};

class ConfinementConstantConcentration : public SOSEvent
{
public:

    ConfinementConstantConcentration(const SOSSolver &solver);

    double newConcentration() const;

private:

    double m_V0;
    double m_c0;

    double m_deltaSum;

    double m_currentVolume;

    // Event interface
public:
    void execute();
    void initialize();
    void reset();
};
