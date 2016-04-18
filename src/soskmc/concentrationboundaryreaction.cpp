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
    m_location(orientation == 0 ? 0 : (dim == 0 ? solver.length() : solver.width())),
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

double ConcentrationBoundaryReaction::freeBoundaryArea() const
{
    const double &hl = solver().confiningSurfaceEvent().height();

    double dhSum = 0;
    for (uint xi = 0; xi < span(); ++xi)
    {
        const int h = heightAtBoundary(xi);

        double dh = hl - h - 2;

        if (dh < 0)
        {
            dh = 0;
        }

        dhSum += dh;
    }

    return dhSum;
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

bool ConcentrationBoundaryReaction::pointIsOnBoundary(const uint x, const uint y) const
{
    if (dim() == 0)
    {
        return x == location();
    }

    else
    {
        return y == location();
    }
}

double ConcentrationBoundaryReaction::_rateExpression(const double freeArea) const
{
    BADAssBool(!solver().concentrationIsZero());

    return m_c*freeArea*solver().diffusionEvent().DScaled();
}

bool ConcentrationBoundaryReaction::isAllowed() const
{
    return freeBoundaryArea() != 0;
}

void ConcentrationBoundaryReaction::executeAndUpdate()
{
    double dhSum = 0;
    const double &hl = solver().confiningSurfaceEvent().height();
    for (uint xi = 0; xi < span(); ++xi)
    {
        const int h = heightAtBoundary(xi);

        double dh = hl - h - 2;

        if (dh < 0)
        {
            dh = 0;
        }

        m_boundaryHeights(xi) = dh;

        dhSum += dh;
        m_accuBoundaryHeights(xi) = dhSum;
    }

    BADAss(dhSum, >=, 0);

    const double R = rng.uniform()*dhSum;

    const uint x0 = binarySearchAndScan(m_accuBoundaryHeights.memptr(),
                                        span(),
                                        R);

    const int h0 = heightAtBoundary(x0);

    const double z0 = h0 + 1 + m_accuBoundaryHeights(x0) - R;

    if (dim() == 0)
    {
        solver().diffusionEvent().executeConcentrationBoundaryReaction(location(), x0, z0);
    }
    else
    {
        solver().diffusionEvent().executeConcentrationBoundaryReaction(x0, location(), z0);
    }
}

double ConcentrationBoundaryReaction::rateExpression()
{
    return _rateExpression(freeBoundaryArea());
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
            const uint &x1 = solver().currentSurfaceChange().x1;
            const uint &y1 = solver().currentSurfaceChange().y1;

            onBoundary = onBoundary || pointIsOnBoundary(x1, y1);
        }

        if (onBoundary)
        {
            m_solver.registerAffectedReaction(this);
        }
    }

    else
    {
        const double heightChange = solver().confiningSurfaceEvent().height() -
                solver().confiningSurfaceEvent().currentConfinementChange().prevHeight;
        changeRate(rate() + _rateExpression(heightChange*span()));
    }
}
