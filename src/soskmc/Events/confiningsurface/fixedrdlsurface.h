#ifndef FIXEDRDLSURFACE_H
#define FIXEDRDLSURFACE_H

#include "fixedsurface.h"
#include "rdlsurface.h"

class FixedRDLSurface : public FixedSurface, public RDLSurface
{
public:
    FixedRDLSurface(SolidOnSolidSolver &solver,
                    const double E0,
                    const double s0,
                    const double ld,
                    const double height);

    ~FixedRDLSurface();

    // Event interface
public:
    void execute()
    {
        return FixedRDLSurface::execute();
    }
    void reset()
    {
        return RDLSurface::reset();
    }

    // ConfiningSurface interface
public:
    void setupInitialConditions()
    {
        RDLSurface::recalculateAllRDLEnergies();
    }

    void registerHeightChange(const uint x, const uint y, std::vector<DiffusionDeposition *> affectedReactions, const uint n)
    {
        (void) affectedReactions;
        (void) n;

        RDLSurface::recalculateRDLEnergy(x, y);

        FixedSurface::registerHeightChange(x, y, affectedReactions, n);
    }

    double confinementEnergy(const uint x, const uint y)
    {
        return RDLSurface::confinementEnergy(x, y);
    }
};

#endif // FIXEDRDLSURFACE_H
