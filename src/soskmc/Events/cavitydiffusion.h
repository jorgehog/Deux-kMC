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

    double localSurfaceSupersaturation(const uint x, const uint y);

    double volume() const;

    const uint &nParticles() const
    {
        return m_particlePositions.n_cols;
    }

    bool isBlockedPosition(const double x, const double y, const double z) const;

    void _dump(const uint frameNumber) const;

private:

    const double m_D0;
    const double m_D;

    const double m_dt0;
    const double m_dt;

    mat m_particlePositions;
    mat m_F;

    void diffuse(const double dt);

};