#include "solidonsolidsolver.h"
#include "solidonsolidreaction.h"
#include "Events/confiningsurface/confiningsurface.h"
#include "Events/cavitydiffusion.h"
#include "../kmcsolver/boundary/boundary.h"


SolidOnSolidSolver::SolidOnSolidSolver(const uint length,
                                       const uint width,
                                       const Boundary *xBoundary,
                                       const Boundary *yBoundary,
                                       const double alpha,
                                       const double mu) :
    KMCSolver({xBoundary, yBoundary}),
    m_dim((( length == 1 ) || ( width == 1 ) ) ? 1 : 2),
    m_length(length),
    m_width(width),
    m_alpha(alpha),
    m_mu(mu),
    m_heights(length, width, fill::zeros),
    m_nNeighbors(length, width),
    m_siteReactions(length, width),
    m_affectedReactions(5)

{
    for (uint x = 0; x < length; ++x)
    {
        for (uint y = 0; y < width; ++y)
        {
            m_siteReactions(x, y) = new DiffusionDeposition(x, y, *this);
        }
    }
}

SolidOnSolidSolver::~SolidOnSolidSolver()
{
    for (uint x = 0; x < length(); ++x)
    {
        for (uint y = 0; y < width(); ++y)
        {
            delete m_siteReactions(x, y);
        }
    }

    m_siteReactions.clear();
}

void SolidOnSolidSolver::registerHeightChange(const uint x, const uint y, const int value)
{
    m_heights(x, y) += value;

    const uint left = leftSite(x);
    const uint right = rightSite(x);
    const uint bottom = bottomSite(y);
    const uint top = topSite(y);

    uint n = 0;

    setNNeighbors(x, y);
    m_affectedReactions.at(n) = &reaction(x, y);
    n++;

    if (!boundary(0)->isBlocked(left))
    {
        setNNeighbors(left, y);
        m_affectedReactions.at(n) = &reaction(left, y);
        n++;
    }

    if (!boundary(0)->isBlocked(right))
    {
        setNNeighbors(right, y);
        m_affectedReactions.at(n) = &reaction(right, y);
        n++;
    }

    if (!boundary(1)->isBlocked(top))
    {
        setNNeighbors(x, top);
        m_affectedReactions.at(n) = &reaction(x, top);
        n++;
    }

    if (!boundary(1)->isBlocked(bottom))
    {
        setNNeighbors(x, bottom);
        m_affectedReactions.at(n) = &reaction(x, bottom);
        n++;
    }

    m_confiningSurfaceEvent->registerHeightChange(x, y, m_affectedReactions, n);

    m_diffusionEvent->registerHeightChange(x, y, value);

    for (uint i = 0; i < n; ++i)
    {
        DiffusionDeposition *reaction = m_affectedReactions.at(i);
        reaction->calculateRate();
    }

}

void SolidOnSolidSolver::setNNeighbors(const uint x, const uint y)
{
    m_nNeighbors(x, y) = calculateNNeighbors(x, y);
}


void SolidOnSolidSolver::setConfiningSurfaceEvent(ConfiningSurface &confiningSurfaceEvent)
{
    m_confiningSurfaceEvent = &confiningSurfaceEvent;
    m_diffusionEvent->setDependency(confiningSurfaceEvent);
}

void SolidOnSolidSolver::setDiffusionEvent(CavityDiffusion &diffusionEvent)
{
    m_diffusionEvent = &diffusionEvent;
}

double SolidOnSolidSolver::confinementEnergy(const uint x, const uint y) const
{
    if (!m_confiningSurfaceEvent->hasStarted())
    {
        return 0;
    }

    return m_confiningSurfaceEvent->confinementEnergy(x, y);
}

double SolidOnSolidSolver::localSurfaceSupersaturation(const uint x, const uint y) const
{
    if (!m_diffusionEvent->hasStarted())
    {
        return 1;
    }

    return m_diffusionEvent->localSurfaceSupersaturation(x, y);
}

uint SolidOnSolidSolver::calculateNNeighbors(const uint x, const uint y) const
{    
    bool connectedLeft = false;
    bool connectedRight = false;
    bool connectedTop = false;
    bool connectedBottom = false;

    uint n = 1;

    const int &h = m_heights(x, y);

    const int left = leftSite(x);
    const int right = rightSite(x);
    const int top = topSite(y);
    const int bottom = bottomSite(y);

    if (!boundary(0)->isBlocked(left))
    {
        const int &hLeft = m_heights(left, y);
        connectedLeft = hLeft >= h;
    }

    if (!boundary(0)->isBlocked(right))
    {
        const int &hRight = m_heights(right, y);
        connectedRight = hRight >= h;
    }

    if (!boundary(1)->isBlocked(top))
    {
        const int &hTop = m_heights(x, top);
        connectedTop = hTop >= h;
    }

    if (!boundary(1)->isBlocked(bottom))
    {
        const int &hBottom = m_heights(x, bottom);
        connectedBottom = hBottom >= h;
    }

    if (connectedLeft)
    {
        n++;
    }

    if (connectedRight)
    {
        n++;
    }

    if (connectedTop)
    {
        n++;
    }

    if (connectedBottom)
    {
        n++;
    }

    return n;

}

int SolidOnSolidSolver::topSite(const uint site, const uint n) const
{
    return boundary(1)->transformCoordinate(site + n);
}

int SolidOnSolidSolver::bottomSite(const uint site, const uint n) const
{
    return boundary(1)->transformCoordinate((int)site - (int)n);
}

int SolidOnSolidSolver::leftSite(const uint site, const uint n) const
{
    return boundary(0)->transformCoordinate((int)site - (int)n);
}

int SolidOnSolidSolver::rightSite(const uint site, const uint n) const
{
    return boundary(0)->transformCoordinate(site + n);
}

uint SolidOnSolidSolver::span() const
{
    int min = m_heights.min();

    if (confiningSurfaceEvent().hasStarted())
    {
        return confiningSurfaceEvent().height() - min;
    }
    else
    {
        return 5*(m_heights.max() - min);
    }
}

void SolidOnSolidSolver::setMu(const double mu)
{
    if (hasStarted())
    {
        double expFac = exp(m_mu - mu);

        for (uint x = 0; x < length(); ++x)
        {
            for (uint y = 0; y < width(); ++y)
            {
                DiffusionDeposition &_reaction = reaction(x, y);
                _reaction.setDiffusionRate(_reaction.dissolutionRate()*expFac);
            }
        }
    }

    m_mu = mu;

}

void SolidOnSolidSolver::initializeSolver()
{
    m_heights.zeros();

    int left;
    int bottom;


    for (uint x = 0; x < m_length; ++x)
    {
        for (uint y = 0; y < m_width; ++y)
        {
            left = leftSite(x);
            bottom = bottomSite(y);

            if (boundary(0)->isBlocked(left) || boundary(1)->isBlocked(bottom))
            {
                continue;
            }

            m_heights(x, y) = (m_heights(left, y) + m_heights(x, bottom))/2 + round(-1 + 2*rng.uniform());
        }
    }

    for (uint x = 0; x < m_length; ++x)
    {
        for (uint y = 0; y < m_width; ++y)
        {
            setNNeighbors(x, y);
        }
    }

    m_confiningSurfaceEvent->setupInitialConditions();
    m_diffusionEvent->setupInitialConditions();

}

uint SolidOnSolidSolver::numberOfReactions() const
{
    return m_siteReactions.size();
}

Reaction *SolidOnSolidSolver::getReaction(const uint n) const
{
    return m_siteReactions(n);
}
