#include <SOSkMC.h>

#include "../apputils.h"

#include "extraneighbor.h"
#include "rdlextraneighborsurface.h"

using ignis::Lattice;

int main()
{
    rng.initialize(time(nullptr));

    const uint L = 100;
    const uint W = 1;

    const double alpha = 2.0;
    const double supersaturation = 0.85;

    const double ld = 4.0;
    const double s0 = 1.0;
    const double Pl = 0.3;

    const double gamma = log(1 + supersaturation);

    SOSSolver solver(L, W, alpha, gamma, true);

    setBoundariesFromIDs(&solver, {0,0,2,2}, L, W);

    RDLPotential rdlpotential(solver, s0, ld);
    solver.addLocalPotential(&rdlpotential);

    ExtraNeighbor extraNeighbor(solver, s0);
    solver.addLocalPotential(&extraNeighbor);

    RDLExtraNeighborSurface rdlSurface(solver, rdlpotential, extraNeighbor, Pl);
    ConstantConcentration constantConcentration(solver);

    const uint interval = 1000;
    DumpSystem dumper(solver, interval);

    Lattice lattice;

    lattice.addEvent(solver);
    lattice.addEvent(rdlSurface);
    lattice.addEvent(constantConcentration);
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

    lattice.eventLoop(1000000);

    cout << rdlSurface.height() << endl;

    return 0;
}
