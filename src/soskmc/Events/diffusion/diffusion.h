#pragma once

#include "../sosevent.h"
#include "../../observers.h"
#include "../../subjects.h"

using kMC::Observer;
using kMC::Subjects;

class SOSDiffusionReaction;
class ConcentrationBoundaryReaction;

class Diffusion : public ignis::LatticeEvent, public Observer<Subjects>
{
public:
    Diffusion(SOSSolver &solver,
              string type,
              string unit = "",
              bool hasOutput = false,
              bool storeValue = false);

    virtual ~Diffusion();

    virtual void dump(const uint frameNumber, const string path = "/tmp", const string ext = "") const;

    virtual double depositionRate(const uint x, const uint y) const = 0;

    virtual bool countPaths() const = 0;

    virtual void executeDiffusionReaction(SOSDiffusionReaction *reaction, const int x, const int y, const int z) = 0;

    virtual void executeConcentrationBoundaryReaction(const uint x, const uint y, const double z) = 0;

    virtual bool isBlockedPosition(const uint x, const uint y, const int z) const = 0;

    virtual double concentration() const = 0;

    virtual bool hasDiscreteParticles() const = 0;

    virtual uint numberOfParticles() const;

    virtual void insertRandomParticle() {}
    virtual void removeRandomParticle() {}

    const SOSSolver &solver() const
    {
        return m_solver;
    }

    double DUnscaled() const;

    double DScaled() const;

private:

    SOSSolver &m_solver;

protected:

    SOSSolver &solver()
    {
        return m_solver;
    }

};
