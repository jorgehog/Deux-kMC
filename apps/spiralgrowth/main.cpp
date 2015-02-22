#include "spiralgrowthsolver.h"

#include <utils.h>

#include <iostream>

#include <time.h>

using namespace ignis;

int main()
{
    rng.initialize(time(NULL));

    const uint L = 100;
    const double alpha = 1.0;

    SolidOnSolidSolver solver(L, alpha);


    Lattice lattice;
    lattice.addEvent(solver);

    lattice.enableOutput(true, 10000);
    lattice.enableProgressReport();
    lattice.eventLoop(1000000);

    return 0;
}

