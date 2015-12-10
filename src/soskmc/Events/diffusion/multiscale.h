#pragma once

#include "offlatticemontecarlo.h"
#include "latticediffusion.h"

class SOSDiffusionReaction;

class Multiscale : public OfflatticeMonteCarlo, public LatticeDiffusion
{
public:
    Multiscale(SOSSolver &solver,
               const double maxdt,
               const uint boundarySpacing = 3);
    ~Multiscale();

    bool checkIfEnoughRoom() const;


private:

    const uint m_boundarySpacing;

    // Event interface
public:
    void execute();
    void reset()
    {
        OfflatticeMonteCarlo::reset();
    }

    // Diffusion interface
public:

    double depositionRate(const uint x, const uint y) const
    {
        return LatticeDiffusion::depositionRate(x, y);
    }

    uint dissolutionPaths(const uint x, const uint y) const
    {
        return OfflatticeMonteCarlo::dissolutionPaths(x, y);
    }

    void executeDiffusionReaction(SOSDiffusionReaction *reaction, const int x, const int y, const int z);

    void executeConcentrationBoundaryReaction(ConcentrationBoundaryReaction *reaction);

    bool isBlockedPosition(const uint x, const uint y, const int z) const
    {
        return LatticeDiffusion::isBlockedPosition(x, y, z);
    }

    void dump(const uint frameNumber, const string path = "/tmp") const;

    double concentration() const;

    bool hasDiscreteParticles() const;
    uint numberOfParticles() const;
    void insertRandomParticle();
    void removeRandomParticle();


    // kMC::Observer interface
public:
    void initializeObserver(const Subjects &subject);
    void notifyObserver(const Subjects &subject)
    {
        LatticeDiffusion::notifyObserver(subject);
    }


    // OfflatticeMonteCarlo interface
public:
    double calculateLocalRateOverD(const uint x, const uint y, const uint n) const;
    double calculateLocalRateOverD(const double rSquared) const {(void)rSquared; return 0;}
};

