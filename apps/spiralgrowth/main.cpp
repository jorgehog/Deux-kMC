#include "solidonsolidsolver.h"

#include "solidonsolidevents.h"

#include <utils.h>

#include <iostream>

#include <time.h>

int main()
{
    rng.initialize(time(NULL));

    const uint L = 100;
    const double alpha = 1.0;

    SolidOnSolidSolver solver(L, alpha);
    SurfaceSize size;
    DumpHeights dumpHeights;

    size.setDependency(solver);
    dumpHeights.setDependency(solver);

    Lattice lattice;
    lattice.addEvent(solver);
    lattice.addEvent(size);
    lattice.addEvent(dumpHeights);

    lattice.enableOutput(true, 10000);
    lattice.enableProgressReport();
    lattice.enableEventValueStorage(true, true, "ignisSOS.ign");
    lattice.eventLoop(1000000);

    return 0;
}

