#include "concentrationboundaryreaction.h"
#include "sossolver.h"
#include "Events/diffusion/diffusion.h"
#include "Events/confiningsurface/confiningsurface.h"


ConcentrationBoundaryReaction::ConcentrationBoundaryReaction(const uint dim,
                                                             const uint orientation,
                                                             SOSSolver &solver,
                                                             const double omega) :
    Reaction(),
    m_dim(dim),
    m_orientation(orientation),
    m_location(orientation == 0 ? 0 : (dim == 0 ? solver.length() - 1 : solver.width() - 1)),
    m_span(dim == 0 ? solver.width() : solver.length()),
    m_solver(solver),
    m_c((omega + 1)*solver.c0()),
    m_boundaryHeights(m_span),
    m_accuBoundaryHeights(m_span)
{
    //location = 0 if orientation = 0 for both dims. If location = 1 we are at the end of the respective dim.
}

ConcentrationBoundaryReaction::~ConcentrationBoundaryReaction()
{

}

uint ConcentrationBoundaryReaction::nBoundarySites() const
{
    const double &hl = solver().confiningSurfaceEvent().height();
    const int hlInt = int(floor(hl));

    return nBoundarySites(hlInt);
}

uint ConcentrationBoundaryReaction::nBoundarySites(const int hlInt) const
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

int ConcentrationBoundaryReaction::heightAtBoundary(const uint n) const
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

bool ConcentrationBoundaryReaction::pointIsOnBoundary(const int x, const int y) const
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

void ConcentrationBoundaryReaction::getBoundaryPosition(uint &xi,
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
}

double ConcentrationBoundaryReaction::_rateExpression(const uint nBoundaryParticles) const
{
    BADAssBool(!solver().concentrationIsZero());

    return m_c*nBoundaryParticles*solver().diffusionEvent().DScaled();
}

void ConcentrationBoundaryReaction::insertParticle(const uint xi, const int h)
{
    if (dim() == 0)
    {
        if (solver().isSurfaceSite(location(), xi, h))
        {
            solver().registerConcentrationBoundaryDeposition(true);
            solver().registerHeightChange(location(), xi, +1);
            solver().registerConcentrationBoundaryDeposition(false);
        }

        else
        {
            solver().diffusionEvent().executeConcentrationBoundaryReaction(location(), xi, h);
        }
    }
    else
    {
        if (solver().isSurfaceSite(xi, location(), h))
        {
            solver().registerConcentrationBoundaryDeposition(true);
            solver().registerHeightChange(xi, location(), +1);
            solver().registerConcentrationBoundaryDeposition(false);
        }

        else
        {
            solver().diffusionEvent().executeConcentrationBoundaryReaction(xi, location(), h);
        }
    }
}

bool ConcentrationBoundaryReaction::isAllowed() const
{
    return nBoundarySites() != 0;
}

void ConcentrationBoundaryReaction::executeAndUpdate()
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

double ConcentrationBoundaryReaction::rateExpression()
{
    return _rateExpression(nBoundarySites());
}


void ConcentrationBoundaryReaction::initializeObserver(const Subjects &subject)
{
    (void) subject;
}

void ConcentrationBoundaryReaction::notifyObserver(const Subjects &subject)
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
