#include <SOSkMC.h>

#include "../apputils.h"

#include "extraneighbor.h"

using ignis::Lattice;

int main()
{
    const uint L = 256;
    const uint W = 1;

    const double alpha = 1.0;
    const double supersaturation = 0.01;
    const double gamma = log(1 + supersaturation);

    SOSSolver solver(L, W, alpha, gamma);

    setBoundariesFromIDs(&solver, {0,0,0,0}, L, W);

    RDLSurface confiningSurface(solver, 0.01*L*W, 0.05, 10.0);
    RDLPotential rdlpotential(solver, confiningSurface);
    solver.addLocalPotential(&rdlpotential);

    ConstantConcentration diffusionEvent(solver);

    solver.setConfiningSurfaceEvent(confiningSurface);
    solver.setDiffusionEvent(diffusionEvent);

    const uint interval = 1000;
    DumpSystem dumper(solver, interval);

    ExtraNeighbor extraNeighbor(solver);

    solver.addLocalPotential(&extraNeighbor);
    solver.registerObserver(&extraNeighbor);

    Lattice lattice;

    lattice.addEvent(solver);
    lattice.addEvent(confiningSurface);
    lattice.addEvent(diffusionEvent);
    lattice.addEvent(dumper);

#ifndef NDEBUG
    RateChecker checker(solver);
    lattice.addEvent(checker);
#endif

    lattice.enableOutput(true, interval);
    lattice.enableProgressReport();
    lattice.enableEventValueStorage(true,
                                    true,
                                    "ignisSOS.ign",
                                    "/tmp",
                                    interval);

    initializeSurface(solver, "random");

    lattice.eventLoop(10000000);

    return 0;
}
