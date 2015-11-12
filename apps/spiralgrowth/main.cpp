#include "solidonsolidsolver.h"
#include "pressurewall.h"

#include "eqmu.h"
#include "equilibriater.h"
#include "miscevents.h"

#include <utils.h>

#include <iostream>

#include <time.h>

using namespace ignis;

string addProcEnding(string filename, string ending, string procending)
{
    stringstream s;
    s << filename << procending << "." << ending;

    return s.str();
}

int main(int argv, char** argc)
{

    string tail;
    string cfgName;

    if (argv != 1)
    {
        int proc = atoi(argc[1]);

        stringstream ending;
        ending << "_" << proc;

        tail = ending.str();

        cfgName = argc[2];
    }
    else
    {
        tail = "";
        cfgName = "infiles/spiralgrowth.cfg";
    }

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

    const uint &shadowingInt = 0;
    const bool shadowing = shadowingInt == 1;

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

    const uint &repeater = getSetting<uint>(root, "repeater");

    SolidOnSolidSolver solver(L, W, alpha, mu, shadowing);
    PressureWall pressureWallEvent(solver, E0, sigma0, r0);

    Time time(solver);

    AverageHeight averageHeight(solver);
    pressureWallEvent.setDependency(averageHeight);

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

    TimeAveragetor nNeighborsAverage(nNeighbors, "avgNeighbors", "", true);
    nNeighborsAverage.setOnsetTime(thermalization);

    GrowthSpeed speed(solver);
    speed.setDependency(averageHeight);
    speed.setOnsetTime(thermalization);

    EqMu eqMu(solver);
    Equilibriater equilibriater(solver, eqMu, nSamplesMuEq, nSamplesMu);

    eqMu.setDependency(nNeighbors);

    MainMesh<uint> lattice;
    lattice.addEvent(solver);
    lattice.addEvent(time);
    lattice.addEvent(averageHeight);
    lattice.addEvent(nNeighbors);
    lattice.addEvent(nNeighborsAverage);
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
     lattice.addEvent(dumpHeightSlice);


#ifndef NDEBUG
    RateChecker checker(solver);
    lattice.addEvent(checker);
#endif

    lattice.enableOutput(ignisOutput, nCyclesPerOutput);
    lattice.enableProgressReport();
    lattice.enableEventValueStorage(storeIgnisData, storeIgnisData, "ignisSOS.ign", path, ignisDataInterval);
    lattice.eventLoop(nCycles);

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

    H5Wrapper::Root h5root(path +  addProcEnding("spiralgrowth", "h5", tail));

    stringstream sizeDesc;
    sizeDesc << L << "x" << W;
    H5Wrapper::Member &sizeRoot = h5root.addMember(sizeDesc.str());

    stringstream potentialDesc;
    potentialDesc <<  "alpha_" << alpha
                   << "_mu_" << solver.mu()
                   << "_E0_" << E0*pressureWall
                   << "_s0_" << sigma0
                   << "_r0_" << r0
                   << "_n_" << repeater;

    H5Wrapper::Member &potentialRoot = sizeRoot.addMember(potentialDesc.str());


    potentialRoot.addData("GrowthSpeed", speed.value());

    potentialRoot.addData("nNeighbors", nNeighbors.value());
    potentialRoot.addData("avgNNeighbors", nNeighborsAverage.value());

    potentialRoot.addData("usewall", pressureWallInt);
    potentialRoot.addData("sigma0", sigma0);
    potentialRoot.addData("r0", r0);
    potentialRoot.addData("E0", E0);
    potentialRoot.addData("rms", rms.value());

    potentialRoot.addData("muEq", muEq);
    potentialRoot.addData("muEqError", muEqError);
    potentialRoot.addData("muShift", muShift);
    potentialRoot.addData("reset", resetInt);
    potentialRoot.addData("useConcEquil", equilibriateInt);
    potentialRoot.addData("size", size.timeAverage());
    potentialRoot.addData("var", var.value());

    potentialRoot.addData("storeIgnisData", storeIgnisDataInt);

    if (storeIgnisData)
    {
        potentialRoot.addData("heightmap", solver.heights());
        potentialRoot.addData("ignisData", lattice.storedEventValues());
        potentialRoot.addData("ignisEventDescriptions", lattice.outputEventDescriptions());
    }

    potentialRoot.addData("randomSeed", seed);

    return 0;
}

