#pragma once

#include "diffusion.h"

namespace Tests
{
class World;
class PathFinder;
}

struct PathFindingJazz
{
    uint xTrans;
    uint yTrans;
    int xEnd;
    int yEnd;
    int zEnd;
};

class OfflatticeMonteCarlo : public virtual Diffusion
{
public:
    OfflatticeMonteCarlo(SOSSolver &solver,
                         const double maxdt,
                         const int depositionBoxHalfSize,
                         string type = "",
                         string unit = "",
                         bool hasOutput = false,
                         bool storeValue = false);

    virtual ~OfflatticeMonteCarlo();

    virtual double calculateLocalRateOverD(const double rSquared) const = 0;

    virtual double calculateLocalRateOverD(const uint x, const uint y, const uint n) const = 0;

    void diffuse(const double dt);

    void diffuseFull(const double dtFull);

    void removeParticle(const uint n);

    void insertParticle(const double x, const double y, const double z);

    void initializeParticleMatrices(const uint N);

    void scan(const uint n, const uint dim, const double dr, const uint maxSteps = 100);

    void scanForDisplacement(const uint n, uint &dim, double &delta, const double stepSize = 0.05);

    const double &maxdt() const
    {
        return m_maxdt;
    }

    const int &depositionBoxHalfSize() const
    {
        return m_depositionBoxHalfSize;
    }

    double acceptanceRatio() const
    {
        return m_accepted/m_trials;
    }

    const long double &trials() const
    {
        return m_trials;
    }

    const long double &accepted() const
    {
        return m_accepted;
    }

    const double &particlePositions(const uint i, const uint j) const
    {
        return m_particlePositions(i, j);
    }

    const mat &particlePositions() const
    {
        return m_particlePositions;
    }

    const uint &nOfflatticeParticles() const
    {
        return m_nParticles;
    }

    void setParticlePosition(const uint i, const uint j, const double value)
    {
        m_particlePositions(i, j) = value;
    }

    const double &localRates(const uint x, const uint y, const uint n) const
    {
        return m_localRates(x, y, n);
    }

    void dumpDiffusingParticles(const uint frameNumber, const string path = "/tmp") const;

    void clearDiffusingParticles();

    void calculateLocalRates();

    double totalParticleDepositionRate(const uint n) const;

    void selectDepositionReactants();
    void selectDepositionReactant(uint &xSelected, uint &ySelected, const uint n);

    bool isInLineOfSight(const uint n, const uint x, const uint y) const;

private:

    const double m_maxdt;
    const int m_depositionBoxHalfSize;
    const int m_boxSize;

    Tests::World *m_world;
    Tests::PathFinder *m_pathFinder;
    vector<PathFindingJazz*> m_pathFindingJazzes;
    uint m_nPathFinds;

    mat m_particlePositions;
    mat m_F;
    cube m_localRates;
    vec m_localRatesForSite;

    uint m_nParticles;

    //tmp
    vec m_selectedDepositionRates;
    umat m_particleDepositionLocations;

    long double m_accepted;
    long double m_trials;

    vec::fixed<6> m_scanDeltas;
    vec::fixed<6> m_scanAbsDeltas;
    vec::fixed<3> m_scanOriginalPositions;

    // Event interface
public:
    void initialize();
    void reset();


    // Diffusion interface
public:
    virtual void dump(const uint frameNumber, const string path = "/tmp") const;
    uint dissolutionPaths(const uint x, const uint y) const;
    virtual double depositionRate(const uint x, const uint y) const;
    virtual void executeDiffusionReaction(SOSDiffusionReaction *reaction, const int x, const int y, const int z);
    virtual void executeConcentrationBoundaryReaction(ConcentrationBoundaryReaction *reaction);
    virtual bool isBlockedPosition(const uint x, const uint y, const int z) const;
    double concentration() const;
    bool hasDiscreteParticles() const;
    uint numberOfParticles() const;
    void insertRandomParticle();
    void removeRandomParticle();

    // Event interface
public:
    void execute();

    // kMC::Observer interface
public:
    void initializeObserver(const Subjects &subject);
    void notifyObserver(const Subjects &subject);
};


