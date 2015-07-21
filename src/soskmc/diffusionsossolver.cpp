#include "diffusionsossolver.h"

DiffusionSOSSolver::~DiffusionSOSSolver()
{
    cout << "derp" << endl;
}

uint DiffusionSOSSolver::numberOfReactions() const
{
    cout << "derp" << endl;
    return SOSSolver::numberOfReactions() + 0;
}

Reaction *DiffusionSOSSolver::getReaction(const uint n) const
{
    if (n < SOSSolver::numberOfReactions())
    {
        return SOSSolver::getReaction(n);
    }

    else
    {
        uint m = n - SOSSolver::numberOfReactions();

        cout << "DERP " << m << endl;
        return NULL;
    }
}

void DiffusionSOSSolver::initializeSolver()
{
    SOSSolver::initializeSolver();

    cout << "DERP" << endl;
}
