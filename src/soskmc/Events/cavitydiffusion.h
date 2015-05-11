#pragma once

#include "solidonsolidevent.h"

class CavityDiffusion : public SolidOnSolidEvent
{
public:

    CavityDiffusion(SolidOnSolidSolver &solver,
                    const double D,
                    const uint pointsPerLatticeUnit = 1);

    ~CavityDiffusion() {}

    void setupInitialConditions();

    void initialize();

    void execute();

    void reset();

    void registerHeightChange(const uint x, const uint y);

    double localSurfaceSupersaturation(const uint x, const uint y);


private:

    const double m_D;
    const uint m_pointsPerLatticeUnit;

    cube m_saturationField;

};
