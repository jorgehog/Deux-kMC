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

    if (argv == 2)
    {
        int proc = atoi(argc[1]);

        stringstream ending;
        ending << "_" << proc;

        tail = ending.str();
    }
    else
    {
        tail = "";
    }

    Config cfg;
    string cfgName = "infiles/" + addProcEnding("spiralgrowth", "cfg", tail);

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

    const uint &shadowingInt = getSetting<uint>(root, "shadowing");
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


    SolidOnSolidSolver solver(L, W, alpha, mu, shadowing);
    PressureWall pressureWallEvent(solver, E0, sigma0, r0);
    AverageHeight averageHeight(solver);
    pressureWallEvent.setDependency(averageHeight);
    pressureWallEvent.setupInitialConditions();

    SurfaceSize size(solver);
    size.setOnsetTime(thermalization);

    DumpHeights3D dumpHeights3D(solver, path);

    NNeighbors nNeighbors(solver);

    EqMu eqMu(solver);
    Equilibriater equilibriater(solver, eqMu, nSamplesMuEq, nSamplesMu);

    eqMu.setDependency(nNeighbors);

    MainMesh<uint> lattice;
    lattice.addEvent(solver);
    lattice.addEvent(averageHeight);
    lattice.addEvent(nNeighbors);

    if (equilibriate)
    {
        lattice.addEvent(eqMu);
        lattice.addEvent(equilibriater);
    }

    if (pressureWall)
    {
        lattice.addEvent(pressureWallEvent);
    }

    //    lattice.addEvent(size);
    //    lattice.addEvent(dumpHeights3D);


#ifndef NDEBUG
    RateChecker checker(solver);
    lattice.addEvent(checker);
#endif

    lattice.enableOutput(true, nCyclesPerOutput);
    lattice.enableProgressReport();
    lattice.enableEventValueStorage(true, true, "ignisSOS.ign", path);
    lattice.eventLoop(nCycles);

    double muEq = 0;
    double muEqError = 0;
    if (equilibriate)
    {
        equilibriater.finalizeAverages();

        muEq = solver.mu();
        muEqError = equilibriater.error();

        if (reset)
        {
            solver.setMu(muEq + muShift);
            lattice.removeEvent(&eqMu);
            lattice.removeEvent(&equilibriater);
            lattice.eventLoop(nCycles);
        }
    }

    H5Wrapper::Root h5root(path +  addProcEnding("spiralgrowth", "h5", tail));



    stringstream sizeDesc;
    sizeDesc << L << "x" << W;
    H5Wrapper::Member &sizeMember = h5root.addMember(sizeDesc.str());

    stringstream potentialDesc;
    potentialDesc <<  "alpha_" << alpha
                  << "_mu_" << solver.mu()
                  << "_E0_" << E0*pressureWall
                  << "_s0_" << sigma0
                  << "_r0_" << r0;

    H5Wrapper::Member &potentialMember = sizeMember.addMember(potentialDesc.str());


    potentialMember.addData("nNeighbors", nNeighbors.value());

    potentialMember.addData("usewall", pressureWallInt);
    potentialMember.addData("sigma0", sigma0);
    potentialMember.addData("r0", r0);
    potentialMember.addData("E0", E0);

    potentialMember.addData("muShift", muShift);
    potentialMember.addData("muEq", muEq);
    potentialMember.addData("muEqError", muEqError);
    potentialMember.addData("shadowing", shadowingInt);
    potentialMember.addData("reset", resetInt);
    potentialMember.addData("useConcEquil", equilibriateInt);
    //    potentialMember.addData("usediffusion", useDiffusionInt);
    //    potentialMember.addData("useisotropicdiffusion", isotropicDiffusionInt);
    //    potentialMember.addData("size", size.value());
    potentialMember.addData("heightmap", solver.heights());
    potentialMember.addData("ignisData", lattice.storedEventValues());
    potentialMember.addData("ignisEventDescriptions", lattice.outputEventDescriptions());

    potentialMember.addData("randomSeed", seed);

    return 0;
}

