#pragma once

#include "solidonsolidevent.h"


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

