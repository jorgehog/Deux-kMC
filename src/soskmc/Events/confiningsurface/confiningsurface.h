#pragma once

#include "../sosevent.h"
#include "../../observers.h"

#include "../../subjects.h"

using kMC::Subject;
using kMC::Observer;
using kMC::Subjects;

struct CurrentConfinementChange
{
    double prevHeight;
};

class ConfiningSurface : public ignis::LatticeEvent, public Observer<Subjects>, public Subject<Subjects>
{
public:
    ConfiningSurface(SOSSolver &solver,
                     string type,
                     string unit = "",
                     bool hasOutput = false,
                     bool storeValue = false);

    virtual ~ConfiningSurface();

    virtual bool hasSurface() const = 0;

    virtual double confinementEnergy(const uint x, const uint y) = 0;

    virtual bool acceptDiffusionMove(const double x0, const double y0, const double z0,
                                     const double x1, const double y1, const double z1) const = 0;

    virtual double diffusionDrift(const double x, const double y, const double z) const = 0;

    const double &height() const
    {
        BADAssBool(hasSurface());

        return m_height;
    }

    void setHeight(const double height);

    const SOSSolver &solver() const
    {
        return m_solver;
    }

    const CurrentConfinementChange &currentConfinementChange() const
    {
        return *m_ccc;
    }

private:

    CurrentConfinementChange *m_ccc;

    double m_height;

    SOSSolver &m_solver;

protected:

    SOSSolver &mutexSolver() const
    {
        return m_solver;
    }
};


