#include "sossolver.h"
#include "dissolutiondeposition.h"
#include "Events/confiningsurface/confiningsurface.h"
#include "Events/diffusion/constantconcentration.h"
#include "../kmcsolver/boundary/boundary.h"


SOSSolver::SOSSolver(const uint length,
                     const uint width,
                     const double alpha,
                     const double mu,
                     const Boundary *xBoundary,
                     const Boundary *yBoundary) :
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
            m_siteReactions(x, y) = new DissolutionDeposition(x, y, *this);
            addReaction(m_siteReactions(x, y));
        }
    }
}

SOSSolver::~SOSSolver()
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

void SOSSolver::registerHeightChange(const uint x, const uint y, const int value)
{
    m_heights(x, y) += value;

    const uint left = leftSite(x);
    const uint right = rightSite(x);
    const uint bottom = bottomSite(y);
    const uint top = topSite(y);

    uint n = 0;

    setNNeighbors(x, y);
    m_affectedReactions.at(n) = &surfaceReaction(x, y);
    n++;

    if (!boundary(0)->isBlocked(left))
    {
        setNNeighbors(left, y);
        m_affectedReactions.at(n) = &surfaceReaction(left, y);
        n++;
    }

    if (!boundary(0)->isBlocked(right))
    {
        setNNeighbors(right, y);
        m_affectedReactions.at(n) = &surfaceReaction(right, y);
        n++;
    }

    if (!boundary(1)->isBlocked(top))
    {
        setNNeighbors(x, top);
        m_affectedReactions.at(n) = &surfaceReaction(x, top);
        n++;
    }

    if (!boundary(1)->isBlocked(bottom))
    {
        setNNeighbors(x, bottom);
        m_affectedReactions.at(n) = &surfaceReaction(x, bottom);
        n++;
    }

    m_confiningSurfaceEvent->registerHeightChange(x, y, m_affectedReactions, n);

    m_diffusionEvent->registerHeightChange(x, y, value);

    //recalcuate rates for neighbor reactions.
    for (uint i = 0; i < n; ++i)
    {
        DissolutionDeposition *reaction = m_affectedReactions.at(i);
        reaction->calculateRate();
    }

}

void SOSSolver::setNNeighbors(const uint x, const uint y)
{
    m_nNeighbors(x, y) = calculateNNeighbors(x, y);
}

void SOSSolver::setHeight(const uint x, const uint y, const int value)
{
    //extremely slow implementation

    int dh = value - height(x, y);

    if (dh == 0)
    {
        return;
    }

    int direction = dh/abs(dh);

    for (int i = 0; i < abs(dh); i += 1)
    {
        registerHeightChange(x, y, direction);
    }

}


void SOSSolver::setConfiningSurfaceEvent(ConfiningSurface &confiningSurfaceEvent)
{
    m_confiningSurfaceEvent = &confiningSurfaceEvent;
}

void SOSSolver::setDiffusionEvent(Diffusion &diffusionEvent)
{
    m_diffusionEvent = &diffusionEvent;
}

double SOSSolver::volume() const
{
    //sum (h_l - h_i) - 1
    return (m_confiningSurfaceEvent->height() - 1)*area() - arma::accu(m_heights);
}

double SOSSolver::confinementEnergy(const uint x, const uint y) const
{
    if (!m_confiningSurfaceEvent->hasStarted())
    {
        return 0;
    }

    return m_confiningSurfaceEvent->confinementEnergy(x, y);
}

double SOSSolver::depositionRate(const uint x, const uint y) const
{
    if (!m_diffusionEvent->hasStarted())
    {
        return 1;
    }

    return m_diffusionEvent->depositionRate(x, y);
}

uint SOSSolver::calculateNNeighbors(const uint x, const uint y) const
{    
    bool connectedLeft, connectedRight,
            connectedTop, connectedBottom;

    findConnections(x, y, connectedLeft, connectedRight,
                    connectedBottom, connectedTop);

    uint n = 1;

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

uint SOSSolver::numberOfSurroundingSolutionSites(const uint x, const uint y) const
{
    return 6u - nNeighbors(x, y);
}


void SOSSolver::getSolutionSite(const uint x, const uint y,
                                int &dx, int &dy, int &dz,
                                const uint siteNumber) const
{
    /*
     *   4
     * 1 0 2
     *   3
     */

    BADAss(siteNumber, <=, numberOfSurroundingSolutionSites(x, y));

    //#0 site is above the (x, y, h+1) site.
    if (siteNumber == 0)
    {
        dx = 0;
        dy = 0;
        dz = 2;

        return;
    }

    //all other sites are besides (x,y, h+1) and has zs = h+1
    else
    {
        dz = 1;
    }

    uint n = 1;

    bool connectedLeft, connectedRight,
            connectedBottom, connectedTop;

    findConnections(x, y, connectedLeft, connectedRight, connectedBottom, connectedTop);

    if (!connectedLeft)
    {
        if (n == siteNumber)
        {
            dx = -1;
            dy =  0;
            return;
        }

        n++;
    }

    if (!connectedRight)
    {
        if (n == siteNumber)
        {
            dx = 1;
            dy = 0;
            return;
        }

        n++;
    }

    if (!connectedBottom)
    {
        if (n == siteNumber)
        {
            dx =  0;
            dy = -1;
            return;
        }

        n++;
    }

    if (!connectedTop)
    {
        if (n == siteNumber)
        {
            dx = 0;
            dy = 1;
            return;
        }
    }

}

int SOSSolver::topSite(const uint site, const uint n) const
{
    return boundary(1)->transformCoordinate(site + n);
}

int SOSSolver::bottomSite(const uint site, const uint n) const
{
    return boundary(1)->transformCoordinate((int)site - (int)n);
}

int SOSSolver::leftSite(const uint site, const uint n) const
{
    return boundary(0)->transformCoordinate((int)site - (int)n);
}

int SOSSolver::rightSite(const uint site, const uint n) const
{
    return boundary(0)->transformCoordinate(site + n);
}

void SOSSolver::findConnections(const uint x,
                                const uint y,
                                const int h,
                                bool &connectedLeft,
                                bool &connectedRight,
                                bool &connectedBottom,
                                bool &connectedTop) const
{
    connectedLeft = false;
    connectedRight = false;
    connectedBottom = false;
    connectedTop = false;

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
}

uint SOSSolver::span() const
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

bool SOSSolver::isBlockedPosition(const double x, const double y, const double z) const
{
    //DERP: Check periodicity

    bool isOutSideBox_x = (x < 0) || (x > length());

    bool isOutSideBox_y = (y < 0) || (y > width());

    bool isOutSideBox_z = (z > confiningSurfaceEvent().height() - 0.5);

    if (isOutSideBox_x || isOutSideBox_y || isOutSideBox_z)
    {
        return true;
    }

    uint X = uint(x);
    uint Y = uint(y);

    if (X == length())
    {
        X = length() - 1;
    }

    if (Y == width())
    {
        Y = width() - 1;
    }

    //DERP: collide in x-y plane
    return z < height(X, Y) + 0.5;
}

bool SOSSolver::isSurfaceSite(const uint x, const uint y, const int z) const
{
    return height(x, y) == z - 1;
}

void SOSSolver::setMu(const double mu)
{
    if (hasStarted())
    {
        double expFac = exp(m_mu - mu);

        for (uint x = 0; x < length(); ++x)
        {
            for (uint y = 0; y < width(); ++y)
            {
                DissolutionDeposition &_reaction = surfaceReaction(x, y);
                _reaction.setDiffusionRate(_reaction.dissolutionRate()*expFac);
            }
        }
    }

    m_mu = mu;

}

void SOSSolver::initializeSolver()
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
