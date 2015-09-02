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

    virtual double calculateLocalRate(const uint x, const uint y, const uint n) const = 0;

    void diffuse(const double dt);

    void removeParticle(const uint n);

    void insertParticle(const double x, const double y, const double z);

    void initializeParticleMatrices(const uint nOfflatticeParticles, const double zMin);

    const double &maxdt() const
    {
        return m_maxdt;
    }

    static constexpr double D()
    {
        return 1.0;
    }

    double acceptanceRatio() const
    {
        return m_accepted/m_trials;
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
        return m_particlePositions.n_cols;
    }

    void setParticlePosition(const uint i, const uint j, const double value)
    {
        m_particlePositions(i, j) = value;
    }

    const double &localRates(const uint x, const uint y, const uint n)
    {
        return m_localRates(x, y, n);
    }

    void dumpDiffusingParticles(const uint frameNumber) const;

protected:

    double &particlePositions(const uint i, const uint j)
    {
        return m_particlePositions(i, j);
    }

private:
    const double m_maxdt;

    mat m_particlePositions;
    mat m_F;
    cube m_localRates;

    long double m_accepted;
    long double m_trials;


    // Event interface
public:
    void initialize();
    void reset();


    // Diffusion interface
public:
    void setupInitialConditions();
    virtual void dump(const uint frameNumber) const;
    uint dissolutionPaths(const uint x, const uint y) const;
    virtual void registerHeightChange(const uint x, const uint y, const int delta);
    virtual double depositionRate(const uint x, const uint y) const;
    virtual void executeDiffusionReaction(SOSDiffusionReaction *reaction, const int x, const int y, const int z);
    virtual void executeConcentrationBoundaryReaction(ConcentrationBoundaryReaction *reaction);
    virtual bool isBlockedPosition(const uint x, const uint y, const int z) const;
};


