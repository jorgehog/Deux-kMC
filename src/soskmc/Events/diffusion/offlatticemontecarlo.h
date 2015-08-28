#pragma once

#include "diffusion.h"

class OfflatticeMonteCarlo : public virtual Diffusion
{
public:
    OfflatticeMonteCarlo(SOSSolver &solver,
                         const double dt,
                         string type = "",
                         string unit = "",
                         bool hasOutput = false,
                         bool storeValue = false);

    virtual ~OfflatticeMonteCarlo();

    virtual void onInsertParticle(const double x, const double y, const double z) = 0;
    virtual void onRemoveParticle(const uint n) = 0;

    void diffuse(const double dt);

    void removeParticle(const uint n);

    void insertParticle(const double x, const double y, const double z);

    void initializeParticleMatrices(const uint nOfflatticeParticles, const double zMin);

    const double &dt() const
    {
        return m_dt;
    }

    const double &D() const
    {
        return m_D;
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

protected:

    double &particlePositions(const uint i, const uint j)
    {
        return m_particlePositions(i, j);
    }

private:

    const double m_D;
    const double m_dt;

    mat m_particlePositions;
    mat m_F;

    long double m_accepted;
    long double m_trials;


    // Event interface
public:
    void initialize();


    // Diffusion interface
public:
    virtual void dump(const uint frameNumber) const;
    uint dissolutionPaths(const uint x, const uint y) const;

};

