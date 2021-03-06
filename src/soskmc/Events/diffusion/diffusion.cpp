#include "diffusion.h"

#include "../confiningsurface/confiningsurface.h"

#include "sossolver.h"

#include "dissolutiondeposition.h"


Diffusion::Diffusion(SOSSolver &solver,
                     string type, string unit, bool hasOutput, bool storeValue) :
    ignis::LatticeEvent(type, unit, hasOutput, storeValue),
    m_solver(solver)
{
    solver.setDiffusionEvent(*this);
}

Diffusion::~Diffusion()
{

}

void Diffusion::dump(const uint frameNumber, const string path, const string ext) const
{
    double h;

    if (solver().confiningSurfaceEvent().hasSurface())
    {
        h = solver().confiningSurfaceEvent().height();
    }

    else
    {
        h = solver().heights().max();
    }

    const int zMin = solver().heights().min();
    int cutfactor = 1000;

    lammpswriter surfacewriter(6, "surfaces" + ext, path);
    surfacewriter.setSystemSize(solver().length(), solver().width(), h, 0, 0, zMin);
    surfacewriter.initializeNewFile(frameNumber);

    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            if (solver().confiningSurfaceEvent().hasSurface())
            {
                surfacewriter << 1
                              << x
                              << y
                              << h
                              << 0.
                              << 0;
            }

            int cut = solver().height(x, y) - cutfactor;

            int start = cut > zMin ? cut : zMin;

            for (int zSurface = start; zSurface <= solver().height(x, y) - 1; ++zSurface)
            {
                surfacewriter << 2
                              << x
                              << y
                              << zSurface
                              << 6
                              << 6;
            }

            surfacewriter << 2
                          << x
                          << y
                          << solver().height(x, y);

            if (solver().surfaceReaction(x, y).isAllowed())
            {
                surfacewriter << solver().totalSurfaceEnergy(x, y)
                              << solver().nNeighbors(x, y);

            }

            else
            {
                surfacewriter << 0.
                              << solver().nNeighbors(x, y);
            }

        }
    }

    surfacewriter.finalize();
}

uint Diffusion::numberOfParticles() const
{
    if (hasDiscreteParticles())
    {
        throw std::runtime_error("Number of particles is not correctly overloaded.");
    }

    return 0;
}

double Diffusion::DScaled() const
{
    return DUnscaled()*solver().timeUnit();
}




