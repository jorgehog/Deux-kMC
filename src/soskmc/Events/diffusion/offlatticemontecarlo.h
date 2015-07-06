#ifndef OFFLATTICEMONTECARLO_H
#define OFFLATTICEMONTECARLO_H

#include "diffusion.h"

class OfflatticeMonteCarlo : public Diffusion
{
public:
    OfflatticeMonteCarlo(SolidOnSolidSolver &solver,
                         const double D,
                         const double dt);

    ~OfflatticeMonteCarlo();

    void dump(const uint frameNumber) const;

    void diffuse(const double dt);

    void removeParticle(const uint n);

    void insertParticle(const double x, const double y, const double z);

    double calculateTimeStep(const double initialCondition, bool calculateDissolutionRate = false);

    double depositionRate(const uint x, const uint y, double timeStep) const;

    double calculateLocalProbability(const uint x, const uint y, const uint n, const double timeStep) const;

    double calculateLocalProbability(const uint x, const uint y, const uint n) const
    {
        return calculateLocalProbability(x, y, n, m_currentTimeStep);
    }

    const double &dt() const
    {
        return m_dt;
    }

    double acceptanceRatio() const
    {
        return m_accepted/m_trials;
    }

    const mat &particlePositions() const
    {
        return m_particlePositions;
    }

    const uint &nParticles() const
    {
        return m_particlePositions.n_cols;
    }


private:

    const double m_D0;
    const double m_D;

    const double m_dt0;
    const double m_dt;

    double m_currentTimeStep;

    static constexpr double m_eps = 0.00001;

    mat m_particlePositions;
    mat m_F;

    cube m_localProbabilities;

    long double m_accepted;
    long double m_trials;


    // Event interface
public:
    void initialize();

    void execute();

    void reset();

    // Diffusion interface
public:
    void setupInitialConditions();

    double depositionRate(const uint x, const uint y) const
    {
        return depositionRate(x, y, m_currentTimeStep);
    }

    void registerHeightChange(const uint x, const uint y, const int delta);

};

#endif // OFFLATTICEMONTECARLO_H
