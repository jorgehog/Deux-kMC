#include "fluxboundaryreaction.h"
#include "sossolver.h"
#include "Events/diffusion/diffusion.h"
#include "Events/confiningsurface/confiningsurface.h"


FluxBoundaryReaction::FluxBoundaryReaction(const uint dim,
                                           const uint orientation,
                                           SOSSolver &solver,
                                           const double scaledFlux) :
    Reaction(),
    m_dim(dim),
    m_orientation(orientation),
    m_location(orientation == 0 ? 0 : (dim == 0 ? solver.length() - 1 : solver.width() - 1)),
    m_span(dim == 0 ? solver.width() : solver.length()),
    m_solver(solver),
    m_fluxOverD(scaledFlux*solver.c0()),
    m_boundaryHeights(m_span),
    m_accuBoundaryHeights(m_span)
{
    //location = 0 if orientation = 0 for both dims. If location = 1 we are at the end of the respective dim.
}

FluxBoundaryReaction::~FluxBoundaryReaction()
{

}

uint FluxBoundaryReaction::nBoundarySites() const
{
    const double &hl = solver().confiningSurfaceEvent().height();
    const int hlInt = int(floor(hl));

    return nBoundarySites(hlInt);
}

uint FluxBoundaryReaction::nBoundarySites(const int hlInt) const
{
    uint nSites = 0;
    for (uint xi = 0; xi < span(); ++xi)
    {
        const int h = heightAtBoundary(xi);

        const uint n = hlInt - h - 1;

        BADAss(h, <, hlInt);

        nSites += n;
    }

    return nSites;
}

int FluxBoundaryReaction::heightAtBoundary(const uint n) const
{
    if (dim() == 0)
    {
        return solver().height(location(), n);
    }

    else
    {
        return solver().height(n, location());
    }
}

bool FluxBoundaryReaction::pointIsOnBoundary(const int x, const int y) const
{
    if (dim() == 0)
    {
        return x == (int)location();
    }

    else
    {
        return y == (int)location();
    }
}

void FluxBoundaryReaction::getBoundaryPosition(uint &xi,
                                               int &h,
                                               const uint whichSite,
                                               const int hlInt) const
{
    uint scan = 0;
    for (uint _xi = 0; _xi < span(); ++_xi)
    {
        for (h = heightAtBoundary(_xi) + 1; h < hlInt; ++h)
        {
            if (scan == whichSite)
            {
                xi = _xi;
                return;
            }

            scan += 1;
        }
    }

    BADAssBreak("unable to find site.");
}

double FluxBoundaryReaction::_rateExpression(const uint boundaryArea) const
{
    BADAssBool(!solver().concentrationIsZero());

    return m_fluxOverD*boundaryArea*solver().diffusionEvent().DScaled();
}

void FluxBoundaryReaction::insertParticle(const uint xi, const int h)
{
    if (dim() == 0)
    {
        if (solver().isSurfaceSite(location(), xi, h))
        {
            solver().registerFluxBoundaryDeposition(true);
            solver().registerHeightChange(location(), xi, +1);
            solver().registerFluxBoundaryDeposition(false);
        }

        else
        {
            solver().diffusionEvent().executeFluxBoundaryReaction(location(), xi, h);
        }
    }
    else
    {
        if (solver().isSurfaceSite(xi, location(), h))
        {
            solver().registerFluxBoundaryDeposition(true);
            solver().registerHeightChange(xi, location(), +1);
            solver().registerFluxBoundaryDeposition(false);
        }

        else
        {
            solver().diffusionEvent().executeFluxBoundaryReaction(xi, location(), h);
        }
    }
}

bool FluxBoundaryReaction::isAllowed() const
{
    return nBoundarySites() != 0;
}

void FluxBoundaryReaction::executeAndUpdate()
{
    const double &hl = solver().confiningSurfaceEvent().height();
    const int hlInt = int(floor(hl));

    const uint nSites = nBoundarySites(hlInt);

    const uint whichSite = rng.uniform()*nSites;

    uint xi;
    int h;
    getBoundaryPosition(xi, h, whichSite, hlInt);

    insertParticle(xi, h);
}

double FluxBoundaryReaction::rateExpression()
{
    return _rateExpression(nBoundarySites());
}


void FluxBoundaryReaction::initializeObserver(const Subjects &subject)
{
    (void) subject;
}

void FluxBoundaryReaction::notifyObserver(const Subjects &subject)
{
    if (subject == Subjects::SOLVER)
    {
        const uint &x = solver().currentSurfaceChange().x;
        const uint &y = solver().currentSurfaceChange().y;

        bool onBoundary = pointIsOnBoundary(x, y);

        if (solver().currentSurfaceChange().type == ChangeTypes::Double)
        {
            const int &x1 = solver().currentSurfaceChange().x1;
            const int &y1 = solver().currentSurfaceChange().y1;

            onBoundary = onBoundary || pointIsOnBoundary(x1, y1);
        }

        if (onBoundary)
        {
            m_solver.registerAffectedReaction(this);
        }
    }

    else
    {
        if (floor(solver().confiningSurfaceEvent().height()) !=
                floor(solver().confiningSurfaceEvent().currentConfinementChange().prevHeight))
        {
            calculateRate();
        }
    }
}
