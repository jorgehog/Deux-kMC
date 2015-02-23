#include "solidonsolidsolver.h"
#include "pressurewall.h"

#include "eqmu.h"
#include "equilibriater.h"
#include "miscevents.h"

#include <utils.h>

#include <iostream>

#include <time.h>

using namespace ignis;

int main()
{
    int seed = time(NULL);
    rng.initialize(seed);

    const uint nCycles = 1000000;
    const uint thermalization = 10000;
    const uint nCyclesPerOutput = 100;

    const uint L = 256;
    const uint W = 1;

    const double alpha = 1.0;
    const double mu = 0.0;

    const bool shadowing = false;

    const double E0 = -0.1*L;
    const double r0 = 1.0;
    const double sigma0 = 1.0;

    const bool equilibriate = false;

    SolidOnSolidSolver solver(L, W, alpha, mu, shadowing);
    PressureWall pressureWallEvent(solver, E0, sigma0, r0);
    AverageHeight averageHeight(solver);
    pressureWallEvent.setDependency(averageHeight);
    pressureWallEvent.setupInitialConditions();

    SurfaceSize size(solver);
    size.setOnsetTime(thermalization);

    DumpHeights3D dumpHeights3D(solver);


    EqMu eqMu(solver);
    Equilibriater equilibriater(solver, eqMu);

    MainMesh<uint> lattice;
    lattice.addEvent(solver);
    lattice.addEvent(averageHeight);
    lattice.addEvent(pressureWallEvent);
//    lattice.addEvent(size);
    lattice.addEvent(dumpHeights3D);

    if (equilibriate)
    {
        lattice.addEvent(eqMu);
        lattice.addEvent(equilibriater);
    }

#ifndef NDEBUG
    RateChecker checker(solver);
    lattice.addEvent(checker);
#endif

    lattice.enableOutput(true, nCyclesPerOutput);
    lattice.enableProgressReport();
    lattice.enableEventValueStorage(true, true, "ignisSOS.ign");
    lattice.eventLoop(nCycles);

    return 0;
}

