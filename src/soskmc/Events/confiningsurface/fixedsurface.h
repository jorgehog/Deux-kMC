#ifndef FIXEDSURFACE_H
#define FIXEDSURFACE_H

#include "confiningsurface.h"

class FixedSurface : public virtual ConfiningSurface
{
public:
    FixedSurface(SolidOnSolidSolver &solver,
                 const double height);

    virtual ~FixedSurface();

    // Event interface
public:
    void execute()
    {
        //pass
    }

    // ConfiningSurface interface
public:
    void setupInitialConditions()
    {
        //pass
    }

    void registerHeightChange(const uint x, const uint y, std::vector<DiffusionDeposition *> affectedReactions, const uint n)
    {
        (void) affectedReactions;
        (void) n;

        if (solver().height(x, y) > height())
        {
            cout << "DERP: crashing in surface.." << endl;
        }
    }

    double confinementEnergy(const uint x, const uint y)
    {
        (void) x;
        (void) y;

        return 0;
    }
};

#endif // FIXEDSURFACE_H
