#pragma once

#include <kMC>
#include "solidonsolidsolver.h"

class SolidOnSolidEvent : public LatticeEvent
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


class DumpHeights : public SolidOnSolidEvent
{
public:

    DumpHeights() :
        SolidOnSolidEvent("height", "", true, true),
        m_filename("/tmp/heighmap.arma")
    {

    }

protected:

    void execute();

private:

    const string m_filename;

};
