#include "solidonsolidsolver.h"

#include "solidonsolidevents.h"

#include <utils.h>

#include <iostream>

#include <time.h>

using namespace ignis;

int main()
{
    rng.initialize(time(NULL));

    const uint L = 30;
    const uint W = 30;

    const double alpha = 1.0;
    const double mu = 1.0;

    const bool shadowing = false;

    SolidOnSolidSolver solver(L, W, alpha, mu, shadowing);

    SurfaceSize size;
    DumpHeights3D dumpHeights3D;
    AverageHeight averageHeight;

    EqMu eqMu;
    Equilibriater equilibriater(solver, eqMu);

    size.setDependency(solver);
    dumpHeights3D.setDependency(solver);
    averageHeight.setDependency(solver);

    eqMu.setDependency(solver);

    MainMesh<uint> lattice;
    lattice.addEvent(solver);
//    lattice.addEvent(size);
    lattice.addEvent(dumpHeights3D);
    lattice.addEvent(averageHeight);

    lattice.addEvent(eqMu);
    lattice.addEvent(equilibriater);

    lattice.enableOutput(true, 100);
    lattice.enableProgressReport();
    lattice.enableEventValueStorage(true, true, "ignisSOS.ign");
    lattice.eventLoop(1000000);

    return 0;
}

