#include "latticediffusion.h"


class LatticeDiffusionConstantN : public LatticeDiffusion
{
public:
    LatticeDiffusionConstantN(SOSSolver &solver);

private:

    uint m_targetNParticles;
    double m_targetSeparation;

    // HeightObserver interface
public:
    void notifyObserver();
    void initializeObserver();
};

