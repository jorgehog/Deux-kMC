#pragma once

#include "solidonsolidevent.h"


class SurfaceSize : public SolidOnSolidEvent
{
public:

    SurfaceSize(const SolidOnSolidSolver &solver) :
        SolidOnSolidEvent(solver, "SurfaceSize", "l0", true, true),
        m_relativeHeightSums(solver.length(), solver.width())
    {

    }

    void initialize();

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

class SurfaceVariance : public SolidOnSolidEvent
{
public:

    SurfaceVariance(const SolidOnSolidSolver &solver) :
        SolidOnSolidEvent(solver, "SurfaceVariance", "", true)
    {

    }

    void initialize()
    {
        m_s2 = 0;
    }

    void execute();

private:

    double m_s2;
};


class AverageHeight : public SolidOnSolidEvent
{
public:

    AverageHeight(const SolidOnSolidSolver &solver) :
        SolidOnSolidEvent(solver, "AverageHeight", "l0", true, true)
    {
        setValue(getValue());
    }

    void execute();

    double getValue() const;
};


class DumpHeightSlice : public SolidOnSolidEvent
{
public:
    DumpHeightSlice(const SolidOnSolidSolver &solver,
                    const uint slicePosition=0,
                    const uint axis=0,
                    const string path = "/tmp",
                    const uint nCyclesperOutput = 1000) :
        SolidOnSolidEvent(solver, "DumpHeightSlice"),
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

class DumpHeights3D : public SolidOnSolidEvent
{
public:

    DumpHeights3D(const SolidOnSolidSolver &solver,
                  const string path = "/tmp") :
        SolidOnSolidEvent(solver, "DumpHeights3D"),
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

    NNeighbors(const SolidOnSolidSolver &solver) :
        SolidOnSolidEvent(solver, "nNeighbors", "", true, true),
        m_nNeighbors(solver.length(), solver.width())
    {

    }

    double value() const
    {
        return m_sum/(cycle() + 1);
    }

    double getLocalValue() const
    {
        return m_localValue/solver().area();
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

    double m_sum;

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


