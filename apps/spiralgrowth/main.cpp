#include "solidonsolidsolver.h"

#include "solidonsolidevents.h"

#include <utils.h>

#include <iostream>

#include <time.h>

using namespace ignis;

int main()
{
    rng.initialize(time(NULL));

    const uint L = 100;
    const uint W = 1;
    const double alpha = 1.0;

    SolidOnSolidSolver solver(L, W, alpha);
    SurfaceSize size;
    DumpHeights3D dumpHeights3D;
    AverageHeight averageHeight;

    size.setDependency(solver);
    dumpHeights3D.setDependency(solver);
    averageHeight.setDependency(solver);

    Lattice lattice;
    lattice.addEvent(solver);
    lattice.addEvent(size);
    lattice.addEvent(dumpHeights3D);
    lattice.addEvent(averageHeight);

    lattice.enableOutput(true, 100);
    lattice.enableProgressReport();
    lattice.enableEventValueStorage(true, true, "ignisSOS.ign");
    lattice.eventLoop(1000000);

    return 0;
}

