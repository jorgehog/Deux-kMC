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
        seed = time(nullptr) % 1000000;
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
    const double &supersaturation = getSetting<double>(root, "supersaturation");
    const uint &zeroConcentration = getSetting<uint>(root, "zeroConcentration");

    const double gamma = log(1 + supersaturation);

    const string &initialSurfaceType = getSetting<string>(root, "initialSurfaceType");

    const Setting &boundarySettings = getSetting(root, "boundaries");
    const uint rightBoundaryID = boundarySettings[0][0];
    const uint leftBoundaryID = boundarySettings[0][1];
    const uint bottomBoundaryID = boundarySettings[1][0];
    const uint topBoundaryID = boundarySettings[1][1];

    const uint averageHeightDepth = getSetting<uint>(root, "averageHeightDepth");

    const int boundaryHeight = getSetting<int>(root, "boundaryHeight");

    const uint &concentrationBoundary = getSetting<uint>(root, "concentrationBoundary");

    const uint &confinementInt = getSetting<uint>(root, "confinement");

    const double &confiningSurfaceHeight = getSetting<double>(root, "confiningSurfaceHeight");
    const double &E0 = getSetting<double>(root, "E0dA")*L*W;
    const double &lD = getSetting<double>(root, "lD");
    const double &sigma0 = getSetting<double>(root, "sigma0");

    const uint &diffuseInt = getSetting<uint>(root, "diffuse");
    const double &maxdt = getSetting<double>(root, "maxdt");
    const double &depRateConstant = getSetting<double>(root, "c");

    const uint &autoCorrelationInt = getSetting<uint>(root, "autocorrelation");
    const uint &xCorrSpan = getSetting<uint>(root, "xCorrSpan");
    const uint &yCorrSpan = getSetting<uint>(root, "yCorrSpan");

    const uint &equilibriateInt = getSetting<uint>(root, "equilibriate");
    const bool equilibriate = equilibriateInt == 1;
    const uint &resetInt = getSetting<uint>(root, "reset");
    const bool reset = resetInt == 1;

    const uint &nSamplesGammaEq = getSetting<uint>(root, "nSamplesGammaEq");
    const uint &nSamplesGamma = getSetting<uint>(root, "nSamplesGamma");

    const double &gammaShift = getSetting<double>(root, "muShift");

    const uint &nCycles = getSetting<uint>(root, "nCycles");
    const uint &thermalization = getSetting<uint>(root, "thermalization");

    const uint &isolateEvents = getSetting<uint>(root, "isolateEvents");
    const Setting &isolatedEvents = getSetting(root, "isolatedEvents");

    const uint &nCyclesPerOutput = getSetting<uint>(root, "nCyclesPerOutput");
    const uint &storeIgnisDataInt = getSetting<uint>(root, "storeIgnisData");
    const bool storeIgnisData = storeIgnisDataInt == 1;
    const bool storeIgnisFile = getSetting<uint>(root, "storeIgnisFile") == 1;
    const uint &ignisDataInterval = getSetting<uint>(root, "ignisDataInterval");
    const bool ignisOutput = getSetting<uint>(root, "ignisOutput") == 1;
    const bool dumpParticles = getSetting<uint>(root, "dumpParticles") == 1;
    const uint &dumpInterval = getSetting<uint>(root, "dumpInterval");

    //---End default config loading
    //---Start ignis environment and solver creation
    Lattice lattice;
    SOSSolver solver(L, W, alpha, gamma);

    if (zeroConcentration)
    {
        solver.setZeroConcentration();
    }

    ConfiningSurface *confiningSurface;
    Diffusion *diffusion;

    setBoundariesFromIDs(&solver, {rightBoundaryID,
                                   leftBoundaryID,
                                   bottomBoundaryID,
                                   topBoundaryID},
                         L, W, boundaryHeight, averageHeightDepth);

    if (concentrationBoundary == 1)
    {
        solver.addConcentrationBoundary(0, Boundary::orientations::FIRST);
    }

    // Selecting confinement model
    if (confinementInt == 1)
    {
        confiningSurface = new NoConfinement(solver);

        if (diffuseInt != 1 && diffuseInt != 5)
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
        confiningSurface = new FixedRDLSurface(solver, sigma0, lD, confiningSurfaceHeight);
    }
    else if (confinementInt == 4)
    {
        confiningSurface = new RDLSurface(solver, E0, sigma0, lD);
    }
    else if (confinementInt == 5)
    {
        confiningSurface = new ConstantConfinement(solver, confiningSurfaceHeight);
    }
    else
    {
        cout << "Invalid confinement: " << confinementInt << endl;
        return 1;
    }
    //

    //selection diffusion model
    if (diffuseInt == 1)
    {
        diffusion = new ConstantConcentration(solver);
    }
    else if (diffuseInt == 2)
    {
        diffusion = new LatticeDiffusion(solver);
    }
    else if (diffuseInt == 3)
    {
        diffusion = new FixedPointTimeStepping(solver, maxdt);
    }
    else if (diffuseInt == 4)
    {
        diffusion = new FirstPassageContinuum(solver, maxdt, depRateConstant);
    }
    else if (diffuseInt == 5)
    {
        diffusion = new ConcentrationProfile(solver, [&L, &gamma] (const uint x, const uint y)
        {
            (void) y;

            static constexpr double c0 = 3.0;
            static constexpr double c1 = 1.0;

            return c0-(c0 - c1)/(L-1.)*x;
        });
    }
    else if (diffuseInt == 6)
    {
        diffusion = new ConfinedConstantConcentration(solver);
    }
    else if (diffuseInt == 7)
    {
        diffusion = new LatticeDiffusionConstantN(solver);
    }
    else
    {
        cout << "Invalid diffusion: " << diffuseInt << endl;
        return 1;
    }

    diffusion->setDependency(confiningSurface);

    lattice.addEvent(solver);
    lattice.addEvent(confiningSurface);
    lattice.addEvent(diffusion);

    //

    ConcentrationTracker conc(solver);

    AverageHeight averageHeight(solver);

    AutoCorrelationHeight autoCorrelation(solver, xCorrSpan, yCorrSpan);
    autoCorrelation.setOnsetTime(thermalization);

    EqGamma eqGamma(solver);
    Equilibriater equilibriater(solver, eqGamma, nSamplesGammaEq, nSamplesGamma);

    DumpSystem systemDumper(solver, dumpInterval, path);
    if (dumpParticles)
    {
        lattice.addEvent(systemDumper);
    }

#ifndef NDEBUG
    RateChecker checker(solver);
    lattice.addEvent(checker);
    eqGamma.setDependency(checker);
#endif

    lattice.enableOutput(ignisOutput, nCyclesPerOutput);
    lattice.enableProgressReport();
    lattice.enableEventValueStorage(storeIgnisData,
                                    storeIgnisFile,
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

    //    DumpHeights3D dumpHeights3D(solver, path);
    //    DumpHeightSlice dumpHeightSlice(solver, 0, 0, path, nCyclesPerOutput);

    NNeighbors nNeighbors(solver);

    GrowthSpeed speed(solver);
    speed.setDependency(averageHeight);

    lattice.addEvent(conc);
    lattice.addEvent(currentTime);
    lattice.addEvent(averageHeight);
    lattice.addEvent(nNeighbors);
    lattice.addEvent(speed);
    lattice.addEvent(rms);
    lattice.addEvent(size);
    lattice.addEvent(var);

    if (autoCorrelationInt == 1)
    {
        lattice.addEvent(autoCorrelation);
    }

    initializeSurface(solver, initialSurfaceType);

    //---End explicit implementation
    //---Running simulation

    EquilibrationOrganizer eqOrg(lattice, equilibriate, reset, true);
    eqOrg.prepare({&eqGamma, &equilibriater}, {&averageHeight, confiningSurface});

    if (isolateEvents != 0)
    {
        vector<string> eventTypes(isolatedEvents.getLength());

        for (int i = 0; i < isolatedEvents.getLength(); ++i)
        {
            string eventType = isolatedEvents[i];

            eventTypes[i] = eventType;
        }

        EventIsolator<uint> ei(lattice);

        ei.isolate(eventTypes);
    }

    lattice.eventLoop(nCycles);

    equilibriater.finalizeAverages();
    double gammaEq = equilibriater.averageGamma();

    if (equilibriate && reset)
    {
        solver.setGamma(gammaEq + gammaShift);
    }

    eqOrg.reset(nCycles);

    //---End simulation
    //---Start data dump

    H5Wrapper::Root h5root(path +  addProcEnding(argv, argc, "spiralgrowth", "h5"));

    stringstream sizeDesc;
    sizeDesc << L << "x" << W;
    H5Wrapper::Member &sizeRoot = h5root.addMember(sizeDesc.str());

    timeval tv;
    gettimeofday(&tv, nullptr);

    __int64_t run_ID = 1000*tv.tv_sec + tv.tv_usec/1000 + 10000000000000u*getProc(argv, argc);

    H5Wrapper::Member &simRoot = sizeRoot.addMember(run_ID);

    simRoot["alpha"] = alpha;
    simRoot["gamma"] = gamma;
    simRoot["supersaturation"] = supersaturation;

    simRoot["initialSurfaceType"] = initialSurfaceType;

    simRoot["rightBoundaryID"] = rightBoundaryID;
    simRoot["leftBoundaryID"] = leftBoundaryID;
    simRoot["bottomBoundaryID"] = bottomBoundaryID;
    simRoot["topBoundaryID"] = topBoundaryID;

    simRoot["averageHeightDepth"] = averageHeightDepth;

    simRoot["boundaryHeight"] = boundaryHeight;

    simRoot["concentrationBoundary"] = concentrationBoundary;

    simRoot["confinement"] = confinementInt;
    simRoot["confiningSurfaceHeight"] = confiningSurfaceHeight;
    simRoot["sigma0"] = sigma0;
    simRoot["r0"] = lD;
    simRoot["E0"] = E0;

    simRoot["diffuse"] = diffuseInt;
    simRoot["maxdt"] = maxdt;
    simRoot["depRateConstant"] = depRateConstant;

    simRoot["autoCorrelationInt"] = autoCorrelationInt;
    simRoot["xCorrSpan"] = xCorrSpan;
    simRoot["yCorrSpan"] = yCorrSpan;

    simRoot["useConcEquil"] = equilibriateInt;
    simRoot["reset"] = resetInt;
    simRoot["gammaEq"] = gammaEq;
    simRoot["gammaEqError"] = equilibriater.error();
    simRoot["gammaShift"] = gammaShift;

    simRoot["storeIgnisData"] = storeIgnisDataInt;

    if (storeIgnisData)
    {
        simRoot["heightmap"] = solver.heights();
        simRoot["ignisData"] = lattice.storedEventValues();
        simRoot["ignisEventDescriptions"] = lattice.outputEventDescriptions();
    }

    simRoot["randomSeed"] = seed;

    //---Start Explicit dumps

    simRoot["GrowthSpeed"] = speed.value();
    simRoot["nNeighbors"] = nNeighbors.value();
    simRoot["rms"] = rms.value();
    simRoot["size"] = size.timeAverage();
    simRoot["var"] = var.value();

    if (autoCorrelationInt == 1)
    {
        simRoot["RACF"] = autoCorrelation.autoCorrelation();
    }

    //---End data dump

    //derp
    //    for (auto & dimBoundary : boundaries)
    //    {
    //        for (auto & boundary : dimBoundary)
    //        {
    //            delete boundary;
    //        }
    //    }

    delete diffusion;
    delete confiningSurface;

    return 0;
}

