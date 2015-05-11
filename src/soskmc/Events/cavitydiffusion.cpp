#include "cavitydiffusion.h"

CavityDiffusion::CavityDiffusion(SolidOnSolidSolver &solver,
                                 const double D,
                                 const uint pointsPerLatticeUnit) :
    SolidOnSolidEvent(solver, "Diffusion"),
    m_D(D),
    m_pointsPerLatticeUnit(pointsPerLatticeUnit)
{
    solver.setDiffusionEvent(*this);
}

void CavityDiffusion::setupInitialConditions()
{
    const uint &res = m_pointsPerLatticeUnit;

    m_saturationField.set_size(solver().length()*res, solver().width()*res, solver().span()*res);
    m_saturationField.fill(1);

    int min = solver().heights().min();
    for (uint x = 0; x < solver().length()*res; ++x)
    {
        for (uint y = 0; y < solver().width()*res; ++y)
        {
            int h = solver().height(x, y);

            for (int z = (h-min)*res + 2; z < solver().span()*res; ++z)
            {
                m_saturationField(x, y, z) = z;
            }

        }
    }

}

void CavityDiffusion::initialize()
{

}

void CavityDiffusion::execute()
{
    stringstream s;
    s << "/tmp/saturfield" << cycle() << ".arma";
    mat xy = m_saturationField(span::all, span(0, 0), span::all);
    xy.save(s.str());
}

void CavityDiffusion::reset()
{
    const uint &res = m_pointsPerLatticeUnit;

    m_saturationField.set_size(solver().length()*res, solver().width()*res, solver().span()*res);

    int min = solver().heights().min();
    for (uint x = 0; x < solver().length()*res; ++x)
    {
        for (uint y = 0; y < solver().width()*res; ++y)
        {
            int h = solver().height(x, y);
            for (uint z = 0; z < h - min + 2; ++z)
            {
                m_saturationField(x, y, z) = 1;
            }
        }
    }
}

void CavityDiffusion::registerHeightChange(const uint x, const uint y)
{
    (void) x;
    (void) y;
}

double CavityDiffusion::localSurfaceSupersaturation(const uint x, const uint y)
{
    int min = solver().heights().min();
    int h = solver().height(x, y);

    return m_saturationField(x, y, h - min);
}

