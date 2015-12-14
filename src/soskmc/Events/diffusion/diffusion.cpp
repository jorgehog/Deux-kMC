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

void Diffusion::dump(const uint frameNumber, const string path) const
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

    lammpswriter surfacewriter(5, "surfaces", path);
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
                              << 0;
            }

            int cut = solver().height(x, y) - 11;

            int start = cut > zMin ? cut : zMin;

            for (int zSurface = start; zSurface <= solver().height(x, y) - 1; ++zSurface)
            {
                surfacewriter << 2
                              << x
                              << y
                              << zSurface
                              << 0;
            }

            surfacewriter << 2
                          << x
                          << y
                          << solver().height(x, y);

            if (solver().surfaceReaction(x, y).isAllowed())
            {
                surfacewriter << solver().surfaceReaction(x, y).dissolutionRate();
            }

            else
            {
                surfacewriter << 0;
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

double Diffusion::DUnscaled() const
{
    return 1.;
}

double Diffusion::DScaled() const
{
    return DUnscaled()*solver().timeUnit();
}




