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
    rng.initialize(time(NULL));

    const uint nCycles = 1000000;
    const uint thermalization = 10000;
    const uint nCyclesPerOutput = 100;

    const uint L = 30;
    const uint W = 30;

    const double alpha = 1.0;
    const double mu = 1.0;

    const bool shadowing = false;

    const double E0 = -1;
    const double r0 = 1.0;
    const double sigma0 = 1.0;

    PressureWall pressureWallEvent(E0, sigma0, r0);
    SolidOnSolidSolver solver(pressureWallEvent, L, W, alpha, mu, shadowing);

    SurfaceSize size;
    size.setOnsetTime(thermalization);

    DumpHeights3D dumpHeights3D;
    AverageHeight averageHeight;

    pressureWallEvent.setDependency(averageHeight);


    EqMu eqMu;
    Equilibriater equilibriater(solver, eqMu);

    size.setDependency(solver);
    dumpHeights3D.setDependency(solver);
    averageHeight.setDependency(solver);

    eqMu.setDependency(solver);

    MainMesh<uint> lattice;
    lattice.addEvent(solver);
    lattice.addEvent(averageHeight);
    lattice.addEvent(pressureWallEvent);
//    lattice.addEvent(size);
    lattice.addEvent(dumpHeights3D);

    lattice.addEvent(eqMu);
    lattice.addEvent(equilibriater);

    lattice.enableOutput(true, nCyclesPerOutput);
    lattice.enableProgressReport();
    lattice.enableEventValueStorage(true, true, "ignisSOS.ign");
    lattice.eventLoop(nCycles);

    return 0;
}

