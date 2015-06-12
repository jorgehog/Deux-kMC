#pragma once

#include "solidonsolidevent.h"

class CavityDiffusion : public SolidOnSolidEvent
{
public:

    CavityDiffusion(SolidOnSolidSolver &solver,
                    const double D);

    ~CavityDiffusion() {}

    void setupInitialConditions();

    void initialize();

    void execute();

    void reset();

    void registerHeightChange(const uint x, const uint y);

    double localSurfaceSupersaturation(const uint x, const uint y);

    double volume() const;

    const uint &nParticles() const
    {
        return m_particlePositions.n_cols;
    }

    bool isBlockedPosition(const double x, const double y, const double z) const;

    void _dump(const uint frameNumber) const;

private:

    const double m_D;

    mat m_particlePositions;
};
