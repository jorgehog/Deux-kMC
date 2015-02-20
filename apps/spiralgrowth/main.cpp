#include "spiralgrowthsolver.h"

#include <utils.h>

#include <iostream>

#include <time.h>

using namespace ignis;

int main()
{

    rng.initialize(time(NULL));

    SpiralGrowthSolver solver(100);

    Lattice lattice;

    lattice.addEvent(solver);

    lattice.enableProgressReport();
    lattice.eventLoop(100000);

    return 0;
}

