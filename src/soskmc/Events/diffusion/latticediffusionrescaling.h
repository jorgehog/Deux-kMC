#include "latticediffusion.h"


class LatticeDiffusionRescaling : public LatticeDiffusion
{
public:
    LatticeDiffusionRescaling(SOSSolver &solver);

private:

    uint m_targetNParticles;
    double m_targetSeparation;

    // HeightConnecter interface
public:
    void registerHeightChange(const uint x, const uint y, const int value, std::vector<DissolutionDeposition *> &affectedSurfaceReactions, const uint nAffectedSurfaceReactions);
    void setupInitialConditions();
};

