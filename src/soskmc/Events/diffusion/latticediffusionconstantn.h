#include "latticediffusion.h"


class LatticeDiffusionConstantN : public LatticeDiffusion
{
public:
    LatticeDiffusionConstantN(SOSSolver &solver);

private:

    uint m_targetNParticles;

    // kMC::Observer interface
public:
    void notifyObserver(const Subjects &subject);
    void initializeObserver(const Subjects &subject);
};

