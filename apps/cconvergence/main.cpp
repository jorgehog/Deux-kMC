#include <SOSkMC.h>

#include "../apputils.h"


int main()
{
    const uint L = 30;
    const uint W = 30;

    const double alpha = 1.0;
    const double gammaEq = 0;


    SOSSolver solver(L, W, alpha, gammaEq);

    setBoundariesFromIDs(&solver, {0,0,0,0}, L, W);

    const double height = 20.0;
    FixedSurface confiningSurface(solver, height);

    const double maxdt = 0.01;
    const double n = 4.0;
    const double c = 1.0;

    FirstPassageContinuum diffusion(solver, maxdt, n, c);

    const uint nSurfacesPerC = 10;
    const uint nSamplesPerSurface = 100;

    double meanDepositionRate;
    double meanDissolutionRate;

    while (true)
    {

        meanDissolutionRate = 0;
        meanDepositionRate = 0;

        for (uint surface = 0; surface < nSurfacesPerC; ++surface)
        {
            initializeSurface(solver, "random");

            solver.initialize();
            confiningSurface.initialize();
            diffusion.initialize();

            for (uint x = 0; x < L; ++x)
            {
                for (uint y = 0; y < W; ++y)
                {
                    meanDissolutionRate += solver.surfaceReaction(x, y).dissolutionRate();
                }
            }


            for (uint sample = 0; sample < nSamplesPerSurface; ++sample)
            {
                diffusion.clearDiffusingParticles();
                diffusion.initializeObserver(Subjects::SOLVER);

                for (uint x = 0; x < L; ++x)
                {
                    for (uint y = 0; y < W; ++y)
                    {
                        meanDepositionRate += diffusion.depositionRate(x, y);
                    }
                }
            }

            cout << "suface " << surface << endl;
        }

        meanDissolutionRate /= (nSurfacesPerC*solver.area());
        meanDepositionRate /= (nSurfacesPerC*solver.area()*nSamplesPerSurface);

        const double ratio = meanDissolutionRate/meanDepositionRate;

        diffusion.setc(diffusion.c()*ratio);

        cout << diffusion.c() << endl;

    }



    return 0;
}
