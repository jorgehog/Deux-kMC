#include "../apputils.h"

#include <SOSkMC.h>

#include <utils.h>

#include <iostream>

#include <time.h>

using namespace ignis;

int main(int argv, char** argc)
{

    //---Start default config loading

    string cfgName = getCfgName(argv, argc);

    Config cfg;
    cfg.readFile(cfgName.c_str());

    const Setting & root = cfg.getRoot();

    const uint &seedType = getSetting<uint>(root, "seedType");

    int seed;

    if (seedType == 0)
    {
        seed = time(NULL);
    }
    else
    {
        seed = getSetting<int>(root, "specificSeed");
    }

    rng.initialize(seed);

    string path = getSetting<string>(root, "path") + "/";

    const uint &L = getSetting<uint>(root, "L");
    const uint &W = getSetting<uint>(root, "W");

    const double &alpha = getSetting<double>(root, "alpha");
    const double &mu = getSetting<double>(root, "mu");

    const uint &pressureWallInt = getSetting<uint>(root, "pressureWall");
    const bool pressureWall = pressureWallInt == 1;

    const double E0 = getSetting<double>(root, "E0dA")*L*W;
    const double &r0 = getSetting<double>(root, "r0");
    const double &sigma0 = getSetting<double>(root, "sigma0");

    const uint &equilibriateInt = getSetting<uint>(root, "equilibriate");
    const bool equilibriate = equilibriateInt == 1;
    const uint &resetInt = getSetting<uint>(root, "reset");
    const bool reset = resetInt == 1;

    const uint &nSamplesMuEq = getSetting<uint>(root, "nSamplesMuEq");
    const uint &nSamplesMu = getSetting<uint>(root, "nSamplesMu");

    const double &muShift = getSetting<double>(root, "muShift");

    const uint &nCycles = getSetting<uint>(root, "nCycles");
    const uint &thermalization = getSetting<uint>(root, "thermalization");
    const uint &nCyclesPerOutput = getSetting<uint>(root, "nCyclesPerOutput");
    const uint &storeIgnisDataInt = getSetting<uint>(root, "storeIgnisData");
    const bool storeIgnisData = storeIgnisDataInt == 1;
    const uint &ignisDataInterval = getSetting<uint>(root, "ignisDataInterval");
    const bool ignisOutput = getSetting<uint>(root, "ignisOutput") == 1;

    //---End default config loading
    //---Start ignis environment and solver creation

    SolidOnSolidSolver solver(L, W, alpha, mu);
    PressureWall pressureWallEvent(solver, E0, sigma0, r0);

    AverageHeight averageHeight(solver);
    pressureWallEvent.setDependency(averageHeight);

    Lattice lattice;
    lattice.addEvent(solver);

#ifndef NDEBUG
    RateChecker checker(solver);
    lattice.addEvent(checker);
#endif

    lattice.enableOutput(ignisOutput, nCyclesPerOutput);
    lattice.enableProgressReport();
    lattice.enableEventValueStorage(storeIgnisData,
                                    storeIgnisData,
                                    "ignisSOS.ign",
                                    path,
                                    ignisDataInterval);


    //---Start Explicit Implementations

    Time currentTime(solver);

    HeightRMS rms(solver);
    rms.setDependency(averageHeight);

    SurfaceSize size(solver);
    size.setOnsetTime(thermalization);

    SurfaceVariance var(solver);
    var.setOnsetTime(thermalization);

    var.setDependency(size);

    DumpHeights3D dumpHeights3D(solver, path);
    DumpHeightSlice dumpHeightSlice(solver, 0, 0, path, nCyclesPerOutput);

    NNeighbors nNeighbors(solver);

    GrowthSpeed speed(solver);
    speed.setDependency(averageHeight);

    EqMu eqMu(solver);
    Equilibriater equilibriater(solver, eqMu, nSamplesMuEq, nSamplesMu);

    eqMu.setDependency(nNeighbors);


    lattice.addEvent(currentTime);
    lattice.addEvent(averageHeight);
    lattice.addEvent(nNeighbors);
    lattice.addEvent(speed);

    double muEq;
    if (equilibriate)
    {
            lattice.addEvent(eqMu);
            lattice.addEvent(equilibriater);
    }
    else
    {
        lattice.addEvent(rms);
        lattice.addEvent(size);
        lattice.addEvent(var);
    }

    if (pressureWall)
    {
        lattice.addEvent(pressureWallEvent);
    }

    //    lattice.addEvent(dumpHeights3D);
    // lattice.addEvent(dumpHeightSlice);


    //---Running simulation
    lattice.eventLoop(nCycles);

    //---Start post run pre dump

    double muEqError;
    if (equilibriate)
    {

        equilibriater.finalizeAverages();

        muEq = solver.mu();

        muEqError = equilibriater.error();

        if (reset)
        {
            if (muShift != 0)
            {
                solver.setMu(muEq + muShift);
            }

            lattice.removeEvent(&eqMu);
            lattice.removeEvent(&equilibriater);
            lattice.addEvent(rms);
            lattice.addEvent(size);
            lattice.addEvent(var);
            lattice.eventLoop(nCycles);
        }
    }

    //---End post run pre dump
    //---Start data dump

    H5Wrapper::Root h5root(path +  addProcEnding(argv, argc, "spiralgrowth", "h5"));

    stringstream sizeDesc;
    sizeDesc << L << "x" << W;
    H5Wrapper::Member &sizeRoot = h5root.addMember(sizeDesc.str());

    H5Wrapper::Member &simRoot = sizeRoot.addMember(time(NULL));

    simRoot.addData("alpha", alpha);
    simRoot.addData("mu", mu);

    simRoot.addData("usewall", pressureWallInt);
    simRoot.addData("reset", resetInt);
    simRoot.addData("useConcEquil", equilibriateInt);

    simRoot.addData("sigma0", sigma0);
    simRoot.addData("r0", r0);
    simRoot.addData("E0", E0);

    simRoot.addData("storeIgnisData", storeIgnisDataInt);

    if (storeIgnisData)
    {
        simRoot.addData("heightmap", solver.heights());
        simRoot.addData("ignisData", lattice.storedEventValues());
        simRoot.addData("ignisEventDescriptions", lattice.outputEventDescriptions());
    }

    simRoot.addData("randomSeed", seed);

    //---Start Explicit dumps

    simRoot.addData("GrowthSpeed", speed.value());
    simRoot.addData("nNeighbors", nNeighbors.value());

    simRoot.addData("rms", rms.value());
    simRoot.addData("muEq", muEq);
    simRoot.addData("muEqError", muEqError);
    simRoot.addData("muShift", muShift);
    simRoot.addData("size", size.timeAverage());
    simRoot.addData("var", var.value());

    //---End data dump

    return 0;
}

