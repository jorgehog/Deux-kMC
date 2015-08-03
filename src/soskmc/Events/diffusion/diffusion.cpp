#include "diffusion.h"

#include "../confiningsurface/confiningsurface.h"

Diffusion::Diffusion(SOSSolver &solver,
                     string type,
                     string unit,
                     bool hasOutput,
                     bool storeValue) :
    SolidOnSolidEvent(solver, type, unit, hasOutput, storeValue)
{
    solver.setDiffusionEvent(*this);
}

Diffusion::~Diffusion()
{

}

void Diffusion::dump(const uint frameNumber) const
{
    const double &h = solver().confiningSurfaceEvent().height();
    const int zMin = solver().heights().min();

    lammpswriter surfacewriter(5, "surfaces", "/tmp");
    surfacewriter.setSystemSize(solver().length(), solver().width(), h, 0, 0, zMin);
    surfacewriter.initializeNewFile(frameNumber);

    for (uint x = 0; x < solver().length(); ++x)
    {
        for (uint y = 0; y < solver().width(); ++y)
        {
            surfacewriter << 1
                          << x
                          << y
                          << h
                          << 5;

            for (int zSurface = zMin; zSurface <= solver().height(x, y) - 1; ++zSurface)
            {
                surfacewriter << 2
                              << x
                              << y
                              << zSurface
                              << 6;
            }

            surfacewriter << 2
                          << x
                          << y
                          << solver().height(x, y)
                          << solver().nNeighbors(x, y);
        }
    }

    surfacewriter.finalize();
}


