#pragma once

#include "../sosevent.h"
#include "../../heightconnecter.h"

class SOSDiffusionReaction;
class ConcentrationBoundaryReaction;

class Diffusion : public ignis::LatticeEvent, public HeightConnecter
{
public:
    Diffusion(SOSSolver &solver,
              string type,
              string unit = "",
              bool hasOutput = false,
              bool storeValue = false);

    virtual ~Diffusion();

    virtual void dump(const uint frameNumber, const string path = "/tmp") const;

    virtual double depositionRate(const uint x, const uint y) const = 0;

    virtual uint dissolutionPaths(const uint x, const uint y) const = 0;

    virtual void executeDiffusionReaction(SOSDiffusionReaction *reaction, const int x, const int y, const int z) = 0;

    virtual void executeConcentrationBoundaryReaction(ConcentrationBoundaryReaction *reaction) = 0;

    virtual bool isBlockedPosition(const uint x, const uint y, const int z) const = 0;

    virtual double concentration() const = 0;

    const SOSSolver &solver() const
    {
        return m_solver;
    }

private:

    SOSSolver &m_solver;

protected:

    SOSSolver &mutexSolver() const
    {
        return m_solver;
    }

};
