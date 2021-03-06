#pragma once

#include "diffusion.h"

#include <functional>

using std::function;

class ConcentrationProfile : public Diffusion
{
public:

    using funcType = function<double(const uint, const uint)>;

    ConcentrationProfile(SOSSolver &solver, funcType profileFunction);
    virtual ~ConcentrationProfile();

private:

    const funcType m_profileFunction;

    // Event interface
public:
    virtual void execute() {}

    // Diffusion interface
public:
    virtual double depositionRate(const uint x, const uint y) const;
    virtual bool countPaths() const;
    virtual void executeDiffusionReaction(SOSDiffusionReaction *reaction, const int x, const int y, const int z);
    virtual void executeFluxBoundaryReaction(const uint x, const uint y, const double z);
    virtual bool isBlockedPosition(const uint x, const uint y, const int z) const;
    double concentration() const;
    bool hasDiscreteParticles() const;

    // kMC::Observer interface
public:
    virtual void initializeObserver(const Subjects &subject)
    {
        (void) subject;
    }

    virtual void notifyObserver(const Subjects &subject);
};

