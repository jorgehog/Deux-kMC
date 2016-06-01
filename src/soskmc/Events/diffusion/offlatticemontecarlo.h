#pragma once

#include "diffusion.h"


class OfflatticeMonteCarlo : public virtual Diffusion
{
public:
    OfflatticeMonteCarlo(SOSSolver &solver,
                         const double maxdt,
                         string type = "",
                         string unit = "",
                         bool hasOutput = false,
                         bool storeValue = false);

    virtual ~OfflatticeMonteCarlo();

    void diffuse(const double dt);

    void diffuseFull(const double dtFull);

    void diffuseSingleParticle(const uint n, const double dt, const double prefac);

    void diffuseSingleParticle(const uint n, const double dt)
    {
        diffuseSingleParticle(n, dt, sqrt(2*DUnscaled()*dt));
    }

    void removeParticle(const uint n);

    void insertParticle(const double x, const double y, const double z);

    void insertRandomParticles(const uint N);

    void scan(const uint n, const uint dim, const double dr, const uint maxSteps = 100);

    void scanForDisplacement(const uint n, uint &dim, double &delta, const double stepSize = 0.05);

    bool noAvailablePositions() const;

    const double &maxdt() const
    {
        return m_maxdt;
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

    void dumpDiffusingParticles(const uint frameNumber, const string path = "/tmp", const string ext = "") const;

    void clearDiffusingParticles();

    virtual void calculateLocalRatesAndUpdateDepositionRates() = 0;

    double totalParticleDepositionRate(const uint n) const;

    bool isInLineOfSight(const uint n, const uint x, const uint y) const;

    void releaseLockedParticles();

    void resetLocalRates(const uint n);

    void findRandomPosition(double &x0, double &y0, double &z0) const;

    double volume() const;

protected:
    cube m_localRates;
    //create cube of size (2l+1 2l+1 N) which translates more.dont saveas many zeros and
    //transformations are fast.

private:

    vector<uint> m_lockedParticles;

    const double m_maxdt;

    mat m_particlePositions;
    mat m_F;
    vec m_accuRatesForSite;

    uint m_nParticles;

    long double m_accepted;
    long double m_trials;

    vec::fixed<6> m_scanDeltas;
    vec::fixed<6> m_scanAbsDeltas;
    vec::fixed<3> m_scanOriginalPositions;

    vector<uint> m_removalQueue;
    uint m_dumpCounter;

    // Event interface
public:
    void initialize();
    void reset();


    // Diffusion interface
public:
    virtual void dump(const uint frameNumber, const string path = "/tmp", const string ext = "") const;
    bool countPaths() const;
    virtual double depositionRate(const uint x, const uint y) const;
    virtual void executeDiffusionReaction(SOSDiffusionReaction *reaction, const int x, const int y, const int z);
    virtual void executeConcentrationBoundaryReaction(const uint x, const uint y, const double z);
    virtual bool isBlockedPosition(const uint x, const uint y, const int z) const;
    double concentration() const;
    bool hasDiscreteParticles() const;
    uint numberOfParticles() const;
    void insertRandomParticle();
    void removeRandomParticle();

    // kMC::Observer interface
public:
    void initializeObserver(const Subjects &subject);
    void notifyObserver(const Subjects &subject);
};


