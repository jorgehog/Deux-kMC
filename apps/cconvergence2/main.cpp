#include <SOSkMC.h>

#include "../apputils.h"


int main()
{
    const uint L = 20;
    const uint W = 20;

    const double alpha = 1.0;
    const double gammaEq = 0;

    rng.initialize(time(nullptr) % 1000000);

    SOSSolver solver(L, W, alpha, gammaEq);

    setBoundariesFromIDs(&solver, {0,0,0,0}, L, W);

    const double height = 10.0;
    FixedSurface confiningSurface(solver, height);

    const double maxdt = 0.1;
    const double c = 0.05;
    FirstPassageContinuum diffusion(solver, maxdt, 3, c);

    Lattice lattice;
        lattice.enableOutput(false);
//    lattice.enableProgressReport();

    lattice.addEvent(solver);
    lattice.addEvent(confiningSurface);
    lattice.addEvent(diffusion);

    const uint nCyclesPerEquilibrium = 1000;
    const uint nCyclesPerC = 3;
    const uint nIterations = 3;

    BasicExecuteEvent<uint> pauser("pauser", [&] (BasicExecuteEvent<uint> *event)
    {
        if (event->cycle() == 0)
        {
            return;
        }

        else if (event->cycle() % nCyclesPerEquilibrium == 0)
        {
            event->stopLoop();
        }
    });

    const uint nLastSize = 500;
    vec nLast(nLastSize);
    BasicExecuteEvent<uint> nLastConc("nLastC", [&] (BasicExecuteEvent<uint> *event)
    {
        if (event->cycle() < nLastSize)
        {
            nLast(event->cycle()) = diffusion.concentration();
        }

        else
        {
            for (uint n = 0; n < nLastSize - 1; ++n)
            {
                nLast(n) = nLast(n + 1);
            }

            nLast(nLastSize - 1) = diffusion.concentration();
        }
    });

    lattice.addEvent(pauser);
    lattice.addEvent(nLastConc);

    lattice.eventLoop(nCyclesPerEquilibrium*nCyclesPerC*nIterations);

    double c0 = solver.concentration();

    for (uint iteration = 0; iteration < nIterations; ++iteration)
    {
        if (iteration != 0)
        {
            lattice.reConnect();
        }

        double ceq = 0;
        for (uint sample = 0; sample < nCyclesPerC; ++sample)
        {
            const double conc = arma::mean(nLast);
            ceq += conc;

            if (sample != nCyclesPerC - 1)
            {
                lattice.reConnect();
            }

            cout << "sample " << sample << endl;
        }

        ceq /= nCyclesPerC;

        diffusion.setc(diffusion.c()*ceq/c0);

        cout << diffusion.c() << " " << c0 << " " << ceq << endl;
//        sleep(5);

//        c0 = ceq;
    }



    return 0;
}
