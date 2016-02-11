#include <SOSkMC.h>

#include "../apputils.h"

#include "extraneighbor.h"
#include "rdlextraneighborsurface.h"

using ignis::Lattice;

int main()
{
    const uint L = 50;
    const uint W = 1;

    const double alpha = 1.0;
    const double supersaturation = 0.01;

    const double ld = 4.0;
    const double s0 = 1.0;
    const double Pl = 0.1;

    const double gamma = log(1 + supersaturation);

    SOSSolver solver(L, W, alpha, gamma);

    setBoundariesFromIDs(&solver, {0,0,0,0}, L, W);

    RDLPotential rdlpotential(solver, s0, ld);
    solver.addLocalPotential(&rdlpotential);

    ExtraNeighbor extraNeighbor(solver, s0);
    solver.addLocalPotential(&extraNeighbor);

    RDLExtraNeighborSurface rdlSurface(solver, rdlpotential, extraNeighbor, Pl);
    ConstantConcentration constantConcentration(solver);

    const uint interval = 1;
    DumpSystem dumper(solver, interval);

    Lattice lattice;

    lattice.addEvent(solver);
    lattice.addEvent(rdlSurface);
    lattice.addEvent(constantConcentration);
//    lattice.addEvent(dumper);

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

    lattice.eventLoop(10000);

    cout << rdlSurface.height() << endl;

    return 0;
}
