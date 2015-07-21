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
    const double &gamma = getSetting<double>(root, "gamma");

    const Setting &boundaries = getSetting(root, "boundaries");
    const uint xBoundaryID = boundaries[0];
    const uint yBoundaryID = boundaries[1];

    const uint &confinementInt = getSetting<uint>(root, "confinement");

    const double &confiningSurfaceHeight = getSetting<double>(root, "confiningSurfaceHeight");
    const double &E0 = getSetting<double>(root, "E0dA")*L*W;
    const double &lD = getSetting<double>(root, "lD");
    const double &sigma0 = getSetting<double>(root, "sigma0");

    const uint &diffuseInt = getSetting<uint>(root, "diffuse");
    const bool diffuse = diffuseInt == 1;
//    const double &D = getSetting<double>(root, "D");
    const double D = 1./6*exp(2*alpha-gamma);
    const double &dt = getSetting<double>(root, "dt");

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

    ConfiningSurface *confiningSurface;
    Diffusion *diffusion;

    Boundary* xBoundary = getBoundaryFromID(xBoundaryID, L);
    Boundary* yBoundary = getBoundaryFromID(yBoundaryID, W);

    SOSSolver solver(L, W, alpha, gamma, xBoundary, yBoundary);

    AverageHeight averageHeight(solver);

    if (confinementInt == 1)
    {
        confiningSurface = new NoConfinement(solver);

        if (diffuse)
        {
            cout << "Diffusion without confinement is not implemented." << endl;
            return 1;
        }
    }
    else if (confinementInt == 2)
    {
        confiningSurface = new FixedSurface(solver, confiningSurfaceHeight);
    }
    else if (confinementInt == 3)
    {
        confiningSurface = new FixedRDLSurface(solver, E0, sigma0, lD, confiningSurfaceHeight);
    }
    else if (confinementInt == 4)
    {
        confiningSurface = new RDLSurface(solver, E0, sigma0, lD);
        confiningSurface->setDependency(averageHeight);
    }
    else
    {
        cout << "Invalid confinement: " << confinementInt << endl;
        return 1;
    }

    if (diffuse)
    {
        diffusion = new OfflatticeMonteCarloNoBoundary(solver, dt);
        diffusion->setDependency(confiningSurface);
    }

    else
    {
        diffusion = new ConstantConcentration(solver);
    }


    EqMu eqMu(solver);
    Equilibriater equilibriater(solver, eqMu, nSamplesMuEq, nSamplesMu);

    Lattice lattice;
    lattice.addEvent(solver);

#ifndef NDEBUG
    RateChecker checker(solver);
    lattice.addEvent(checker);
    eqMu.setDependency(checker);
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

    lattice.addEvent(currentTime);
    lattice.addEvent(averageHeight);
    lattice.addEvent(nNeighbors);
    lattice.addEvent(speed);
    lattice.addEvent(rms);
    lattice.addEvent(size);
    lattice.addEvent(var);


    lattice.addEvent(confiningSurface);
    lattice.addEvent(diffusion);

    const bool dumpHeights = false;
    if (dumpHeights)
    {
        lattice.addEvent(dumpHeights3D);
        lattice.addEvent(dumpHeightSlice);
    }

    //---End explicit implementation
    //---Running simulation

    EquilibrationOrganizer eqOrg(lattice, equilibriate, reset, true);
    eqOrg.prepare({&eqMu, &equilibriater}, {&averageHeight, confiningSurface});

    lattice.eventLoop(nCycles);

    equilibriater.finalizeAverages();
    double muEq = equilibriater.averageMu();

    solver.setMu(muEq + muShift);

    eqOrg.reset(nCycles);

    //---End simulation
    //---Start data dump

    H5Wrapper::Root h5root(path +  addProcEnding(argv, argc, "spiralgrowth", "h5"));

    stringstream sizeDesc;
    sizeDesc << L << "x" << W;
    H5Wrapper::Member &sizeRoot = h5root.addMember(sizeDesc.str());

    H5Wrapper::Member &simRoot = sizeRoot.addMember(time(NULL));

    simRoot.addData("alpha", alpha);
    simRoot.addData("gamma", gamma);

    simRoot.addData("xBoundaryID", xBoundaryID);
    simRoot.addData("yBoundaryID", yBoundaryID);

    simRoot.addData("confinement", confinementInt);
    simRoot.addData("sigma0", sigma0);
    simRoot.addData("r0", lD);
    simRoot.addData("E0", E0);

    simRoot.addData("diffuse", diffuseInt);
    simRoot.addData("D", D);

    simRoot.addData("useConcEquil", equilibriateInt);
    simRoot.addData("reset", resetInt);
    simRoot.addData("muEq", muEq);
    simRoot.addData("muEqError", equilibriater.error());
    simRoot.addData("muShift", muShift);

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
    simRoot.addData("size", size.timeAverage());
    simRoot.addData("var", var.value());

    //---End data dump

    delete xBoundary;
    delete yBoundary;
    delete diffusion;
    delete confiningSurface;

    return 0;
}

