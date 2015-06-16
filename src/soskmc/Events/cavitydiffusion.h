#pragma once

#include "solidonsolidevent.h"

class CavityDiffusion : public SolidOnSolidEvent
{
public:

    CavityDiffusion(SolidOnSolidSolver &solver,
                    const double D,
                    const double dt);

    ~CavityDiffusion() {}

    void setupInitialConditions();

    void initialize();

    void execute();

    void reset();

    void registerHeightChange(const uint x, const uint y, const int value);

    double localSurfaceSupersaturation(const uint x, const uint y)
    {
        return localSurfaceSupersaturation(x, y, m_currentTimeStep);
    }

    double localSurfaceSupersaturation(const uint x, const uint y, const double timeStep);

    double volume() const;

    const uint &nParticles() const
    {
        return m_particlePositions.n_cols;
    }

    bool isBlockedPosition(const double x, const double y, const double z) const;

    void _dump(const uint frameNumber) const;

    const mat &particlePositions() const
    {
        return m_particlePositions;
    }

    void diffuse(const double dt);

    void removeParticle(const uint n);

    void insertParticle(const double x, const double y, const double z);

    double calculateTimeStep(const double initialCondition, bool calculateDissolutionRate = false);

    const double &dt() const
    {
        return m_dt;
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

};
