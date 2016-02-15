#include <SOSkMC.h>

#include "../apputils.h"

#include "extraneighbor.h"
#include "rdlextraneighborsurface.h"

using ignis::Lattice;

int main()
{
//    rng.initialize(time(nullptr));
    rng.initialize(10000);

    const uint L = 100;
    const uint W = 1;

    const double alpha = 1.0;

    const double ld = 5.0;
    const double s0 = 0.5;
    const double Pl = 0.2;

    const double gamma = alpha*Pl/(ld*(1-exp(-1/ld)));

    SOSSolver solver(L, W, alpha, gamma, true);

    setBoundariesFromIDs(&solver, {0,0,0,0}, L, W);

    RDLPotential rdlpotential(solver, s0, ld);
    solver.addLocalPotential(&rdlpotential);

    ExtraNeighbor extraNeighbor(solver, s0/(L*W));
    solver.addLocalPotential(&extraNeighbor);

    RDLExtraNeighborSurface rdlSurface(solver, rdlpotential, extraNeighbor, Pl);
    ConstantConcentration constantConcentration(solver);

    const uint interval = 1;
    DumpSystem dumper(solver, interval);

    Lattice lattice;

    lattice.addEvent(solver);
    lattice.addEvent(rdlSurface);
    lattice.addEvent(constantConcentration);
    lattice.addEvent(dumper);

    AverageHeight avgH(solver);
    lattice.addEvent(avgH);

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

    rdlSurface.setHeight(rdlSurface.getRdlEquilibrium());

    lattice.eventLoop(1000);

    cout << rdlSurface.height() << endl;

    return 0;
}
