#include "sossolver.h"
#include "dissolutiondeposition.h"
#include "concentrationboundaryreaction.h"
#include "Events/confiningsurface/confiningsurface.h"
#include "Events/diffusion/constantconcentration.h"
#include "../kmcsolver/boundary/boundary.h"
#include "heightconnecter.h"

SOSSolver::SOSSolver(const uint length,
                     const uint width,
                     const double alpha,
                     const double mu) :
    KMCSolver(),
    m_heights_set(false),
    m_dim((( length == 1 ) || ( width == 1 ) ) ? 1 : 2),
    m_confiningSurfaceEvent(nullptr),
    m_diffusionEvent(nullptr),
    m_length(length),
    m_width(width),
    m_alpha(alpha),
    m_mu(mu),
    m_heights(length, width, fill::zeros),
    m_nNeighbors(length, width),
    m_siteReactions(length, width),
    m_affectedSurfaceReactions(5)

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

    for (ConcentrationBoundaryReaction *r : m_concentrationBoundaryReactions)
    {
        delete r;
    }

    m_siteReactions.clear();
    m_concentrationBoundaryReactions.clear();
}

void SOSSolver::registerHeightChange(const uint x, const uint y, const int value)
{
    BADAssBool(!isOutsideBox(x, y));

    m_heights(x, y) += value;

    registerChangedSite(x, y);

    m_averageHeight += value/double(area());

    BADAssClose(averageHeight(), accu(heights())/double(area()), 1E-3);

    const uint left = leftSite(x, y, m_heights(x, y));
    const uint right = rightSite(x, y, m_heights(x, y));
    const uint bottom = bottomSite(x, y, m_heights(x, y));
    const uint top = topSite(x, y, m_heights(x, y));

    uint n = 0;

    m_affectedSurfaceReactions.at(n) = &surfaceReaction(x, y);
    n++;

    if (!isOutsideBoxSingle(left, 0))
    {
        setNNeighbors(left, y);
        m_affectedSurfaceReactions.at(n) = &surfaceReaction(left, y);
        n++;
    }

    if (!isOutsideBoxSingle(right, 0))
    {
        setNNeighbors(right, y);
        m_affectedSurfaceReactions.at(n) = &surfaceReaction(right, y);
        n++;
    }

    if (!isOutsideBoxSingle(top, 1))
    {
        setNNeighbors(x, top);
        m_affectedSurfaceReactions.at(n) = &surfaceReaction(x, top);
        n++;
    }

    if (!isOutsideBoxSingle(bottom, 1))
    {
        setNNeighbors(x, bottom);
        m_affectedSurfaceReactions.at(n) = &surfaceReaction(x, bottom);
        n++;
    }

    for (HeightConnecter *heightConnecter : m_heightConnecters)
    {
        heightConnecter->registerHeightChange(x, y, value, m_affectedSurfaceReactions, n);
    }

    //recalcuate rates for neighbor reactions.
    for (uint i = 1; i < n; ++i)
    {
        registerAffectedReaction(m_affectedSurfaceReactions.at(i));
    }

    updateConcentrationBoundaryIfOnBoundary(x, y);
}

void SOSSolver::registerChangedSite(const uint x, const uint y)
{
    m_changedSurfaceSites.insert(make_pair(x, y));
    setNNeighbors(x, y);
    registerAffectedReaction(&surfaceReaction(x, y));
}

void SOSSolver::setNNeighbors(const uint x, const uint y)
{
    m_nNeighbors(x, y) = calculateNNeighbors(x, y);
}

void SOSSolver::setHeight(const uint x, const uint y, const int value, const bool iteratively)
{
    if (!m_heights_set)
    {
        m_heights_set = true;
    }

    if (!iteratively)
    {
        m_heights(x, y) = value;
        return;
    }

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

    BADAssEqual(m_heights(x, y), value);

}

void SOSSolver::setHeights(const imat &heights, const bool iteratively)
{
    for (uint x = 0; x < length(); ++x)
    {
        for (uint y = 0; y < width(); ++y)
        {
            setHeight(x, y, heights(x, y), iteratively);
        }
    }
}


void SOSSolver::setConfiningSurfaceEvent(ConfiningSurface &confiningSurfaceEvent)
{
    m_confiningSurfaceEvent = &confiningSurfaceEvent;
    registerHeightConnecter(&confiningSurfaceEvent);
}

void SOSSolver::setDiffusionEvent(Diffusion &diffusionEvent)
{
    m_diffusionEvent = &diffusionEvent;
    registerHeightConnecter(&diffusionEvent);
}

double SOSSolver::volume() const
{
    //sum (h_l - h_i - 1)
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

uint SOSSolver::calculateNNeighbors(const uint x, const uint y, const int h) const
{    
    bool connectedLeft, connectedRight,
            connectedTop, connectedBottom;

    findConnections(x, y, h, connectedLeft, connectedRight,
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

uint SOSSolver::numberOfSurroundingSolutionSites(const uint x, const uint y, const int h) const
{
    bool connectedLeft, connectedRight,
            connectedTop, connectedBottom;
    bool connectedLeft2, connectedRight2,
            connectedTop2, connectedBottom2;

    findConnections(x, y, h, connectedLeft, connectedRight,
                    connectedBottom, connectedTop, false);
    findConnections(x, y, h-1, connectedLeft2, connectedRight2,
                    connectedBottom2, connectedTop2);

    bool connectedAbove = isBlockedPosition(x, y, h+1) || diffusionEvent().isBlockedPosition(x, y, h+1);

    uint n = 5;

    if (connectedLeft || connectedLeft2)
    {
        n--;
    }

    if (connectedRight || connectedRight2)
    {
        n--;
    }

    if (connectedBottom || connectedBottom2)
    {
        n--;
    }

    if (connectedTop || connectedTop2)
    {
        n--;
    }

    if (connectedAbove)
    {
        n--;
    }

    return n;

}


void SOSSolver::getSolutionSite(const uint x, const uint y,
                                const int height,
                                int &dx, int &dy, int &dz,
                                const uint siteNumber) const
{
    /*
     *   4
     * 1 0 2
     *   3
     */

    BADAss(siteNumber, <=, numberOfSurroundingSolutionSites(x, y, height));

    bool connectedAbove = isBlockedPosition(x, y, height + 1);

    uint n = 0;

    //#0 site is above the (x, y, h) site.
    if (!connectedAbove)
    {
        if (n == siteNumber)
        {
            dx = 0;
            dy = 0;
            dz = 1;
            return;
        }

        n++;
    }

    bool connectedLeft, connectedRight,
            connectedTop, connectedBottom;
    bool connectedLeft2, connectedRight2,
            connectedTop2, connectedBottom2;

    findConnections(x, y, height, connectedLeft, connectedRight,
                    connectedBottom, connectedTop, false);
    findConnections(x, y, height-1, connectedLeft2, connectedRight2,
                    connectedBottom2, connectedTop2);

    if (!connectedLeft && !connectedLeft2)
    {
        if (n == siteNumber)
        {
            dx = -1;
            dy =  0;
            dz = 0;
            return;
        }

        n++;
    }

    if (!connectedRight && !connectedRight2)
    {
        if (n == siteNumber)
        {
            dx = 1;
            dy = 0;
            dz = 0;
            return;
        }

        n++;
    }

    if (!connectedBottom && !connectedBottom2)
    {
        if (n == siteNumber)
        {
            dx =  0;
            dy = -1;
            dz = 0;
            return;
        }

        n++;
    }

    if (!connectedTop && !connectedTop2)
    {
        if (n == siteNumber)
        {
            dx = 0;
            dy = 1;
            dz = 0;
            return;
        }
    }

}

int SOSSolver::topSite(const uint x, const uint y, const int z, const uint n) const
{
    return boundary(1, 1)->transformCoordinate(y + n, x, z);
}

int SOSSolver::bottomSite(const uint x, const uint y, const int z, const uint n) const
{
    return boundary(1, 0)->transformCoordinate((int)y - (int)n, x, z);
}

int SOSSolver::leftSite(const uint x, const uint y, const int z, const uint n) const
{
    return boundary(0, 0)->transformCoordinate((int)x - (int)n, y, z);
}

int SOSSolver::rightSite(const uint x, const uint y, const int z, const uint n) const
{
    return boundary(0, 1)->transformCoordinate(x + n, y, z);
}

void SOSSolver::findConnections(const uint x,
                                const uint y,
                                const int h,
                                bool &connectedLeft,
                                bool &connectedRight,
                                bool &connectedBottom,
                                bool &connectedTop,
                                bool onlySurface) const
{

    const int left = leftSite(x, y, h);
    const int right = rightSite(x, y, h);
    const int top = topSite(x, y, h);
    const int bottom = bottomSite(x, y, h);

    connectedLeft = findSingleConnection(left, 0, 0, y, h, onlySurface);
    connectedRight = findSingleConnection(right, 0, 1, y, h, onlySurface);
    connectedBottom = findSingleConnection(bottom, 1, 0, x, h, onlySurface);
    connectedTop = findSingleConnection(top, 1, 1, x, h, onlySurface);

    //    connectedLeft = false;
    //    connectedRight = false;
    //    connectedBottom = false;
    //    connectedTop = false;
    //    bool checkConnection;

    //    checkConnection = !isOutsideBoxSingle(left, 0);
    //    if (!boundary(0, 0)->isBlocked(left))
    //    {
    //        if (!onlySurface)
    //        {
    //            if (diffusionEvent().isBlockedPosition(left, y, h))
    //            {
    //                connectedLeft = true;
    //                checkConnection = false;
    //            }
    //        }

    //        if (checkConnection)
    //        {
    //            const int &hLeft = m_heights(left, y);
    //            connectedLeft = hLeft >= h;
    //        }
    //    }

    //    else
    //    {
    //        connectedLeft = true;
    //    }

    //    checkConnection = !isOutsideBoxSingle(right, 0);
    //    if (!boundary(0, 1)->isBlocked(right))
    //    {

    //        if (!onlySurface)
    //        {
    //            if (diffusionEvent().isBlockedPosition(right, y, h))
    //            {
    //                connectedRight = true;
    //                checkConnection = false;
    //            }
    //        }

    //        if (checkConnection)
    //        {
    //            const int &hRight = m_heights(right, y);
    //            connectedRight = hRight >= h;
    //        }
    //    }

    //    else
    //    {
    //        connectedRight = true;
    //    }

    //    checkConnection = !isOutsideBoxSingle(top, 1);
    //    if (!boundary(1, 1)->isBlocked(top))
    //    {
    //        if (!onlySurface)
    //        {
    //            if (diffusionEvent().isBlockedPosition(x, top, h))
    //            {
    //                connectedTop = true;
    //                checkConnection = false;
    //            }
    //        }

    //        if (checkConnection)
    //        {
    //            const int &hTop = m_heights(x, top);
    //            connectedTop = hTop >= h;
    //        }
    //    }

    //    else
    //    {
    //        connectedTop = true;
    //    }

    //    checkConnection = !isOutsideBoxSingle(bottom, 1);
    //    if (!boundary(1, 0)->isBlocked(bottom))
    //    {
    //        if (!onlySurface)
    //        {
    //            if (diffusionEvent().isBlockedPosition(x, bottom, h))
    //            {
    //                connectedBottom = true;
    //                checkConnection = false;
    //            }
    //        }

    //        if (checkConnection)
    //        {
    //            const int &hBottom = m_heights(x, bottom);
    //            connectedBottom = hBottom >= h;
    //        }
    //    }

    //    else
    //    {
    //        connectedBottom = true;
    //    }
}

bool SOSSolver::findSingleConnection(const int xNeighbor,
                                     const uint dim,
                                     const uint orientation,
                                     const uint y,
                                     const int h,
                                     bool onlySurface) const
{

    bool connected = false;
    bool checkConnection = !isOutsideBoxSingle(xNeighbor, dim);

    if (!boundary(dim, orientation)->isBlocked(xNeighbor, y, h))
    {
        if (!onlySurface)
        {
            bool diffBlocked;

            if (dim == 0)
            {
                diffBlocked = diffusionEvent().isBlockedPosition(xNeighbor, y, h);
            }

            else
            {
                diffBlocked = diffusionEvent().isBlockedPosition(y, xNeighbor, h);
            }

            if (diffBlocked)
            {
                connected = true;
                checkConnection = false;
            }
        }

        if (checkConnection)
        {
            int hNeighbor;

            if (dim == 0)
            {
                hNeighbor = m_heights(xNeighbor, y);
            }

            else
            {
                hNeighbor = m_heights(y, xNeighbor);
            }

            connected = hNeighbor >= h;
        }
    }

    else
    {
        connected = true;
    }

    return connected;
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

uint SOSSolver::boundaryOrientation(const double x, const uint dim) const
{
    const uint lx = dim == 0 ? length() : width();

    return x >= lx/2 ? 1 : 0;
}

double SOSSolver::boundaryTransform(const double x, const double y, const double z, const uint dim) const
{
    if (dim == 0)
    {
        return closestBoundary(x, dim)->transformCoordinate(x, y, z);
    }

    else if (dim == 1)
    {
        return closestBoundary(x, dim)->transformCoordinate(y, x, z);
    }

    else
    {
        return z;
    }

}

double SOSSolver::boundaryTransform(const double x, const double y, const double z, const double dxi, const uint dim) const
{

    if (dim == 0)
    {
        return closestBoundary(x, dim)->transformCoordinate(x + dxi, y, z);
    }

    else if (dim == 1)
    {
        return closestBoundary(x, dim)->transformCoordinate(y + dxi, x, z);
    }

    else
    {
        return z + dxi;
    }

}

const Boundary *SOSSolver::closestBoundary(const double x, const uint dim) const
{
    return boundary(dim, boundaryOrientation(x, dim));
}

void SOSSolver::addConcentrationBoundary(const uint dim, const Boundary::orientations orientation)
{
    uint orientationInt;
    if (orientation == Boundary::orientations::FIRST)
    {
        orientationInt = 0;
    }

    else
    {
        orientationInt = 1;
    }

    const Boundary *b = boundary(dim, orientationInt);

    ConcentrationBoundaryReaction *concReaction = new ConcentrationBoundaryReaction(dim, orientationInt, *this);


    int outSide = concReaction->location() + (2*orientationInt - 1);
    int outSideTrans;
    bool blocked;

    if (dim == 0)
    {
        outSideTrans = b->transformCoordinate(outSide, 0, height(outSide, 0));
        blocked = b->isBlocked(outSideTrans, 0, height(outSideTrans, 0));
    }

    else if (dim == 1)
    {
        outSideTrans = b->transformCoordinate(outSide, 0, height(0, outSide));
        blocked = b->isBlocked(outSideTrans, 0, height(0, outSideTrans));
    }

    else
    {
        throw std::logic_error("invalid dimension");
    }


    if (blocked)
    {
        throw std::logic_error("Concentration boundaries must not be blocked.");
    }

    if (outSideTrans != outSide)
    {
        throw std::logic_error("Concentration boundaries must not transform coordinates.");
    }

    m_concentrationBoundaryReactions.push_back(concReaction);

    addReaction(concReaction);
}

bool SOSSolver::isBlockedPosition(const double x, const double y, const double z) const
{

    //check vs length and not length - 1 or - 0.5 because of transformation issues.
    //so i.e. for periodicity,
    bool isOutSideBox_x = (x < 0) || (x >= length());

    bool isOutSideBox_y = (y < 0) || (y >= width());

    bool isOutSideBox_z = (z > confiningSurfaceEvent().height() - 0.5) && confiningSurfaceEvent().hasSurface();

    if (isOutSideBox_x || isOutSideBox_y || isOutSideBox_z)
    {
        return true;
    }

    uint X = uint(round(x));
    uint Y = uint(round(y));

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

bool SOSSolver::isOutsideBoxSingle(const int x, const uint dim) const
{
    if (x < 0)
    {
        return true;
    }

    uint ux = (uint)x;

    if (dim == 0)
    {
        return ux >= length();
    }

    else if (dim == 1)
    {
        return ux >= width();
    }

    else
    {
        return false;
    }
}

bool SOSSolver::isOutsideBox(const int x, const int y) const
{
    return isOutsideBoxSingle(x, 0) || isOutsideBoxSingle(y, 1);
}

bool SOSSolver::isSurfaceSite(const uint x, const uint y, const int z) const
{
    if (isOutsideBox(x, y))
    {
        return false;
    }

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

void SOSSolver::updateConcentrationBoundaryIfOnBoundary(const uint x, const uint y)
{
    for (ConcentrationBoundaryReaction *r : m_concentrationBoundaryReactions)
    {
        if (r->pointIsOnBoundary(x, y))
        {
            r->calculateRate();
        }
    }
}

double SOSSolver::closestSquareDistance(const uint x, const uint y, const int z,
                                        const double xp, const double yp, const double zp) const
{

    double dxi0, dxj0, dxk0;
    double dxi1, dxj1, dxk1;

    closestBoundary(x, 0)->closestImage(x, y, z, xp, yp, zp, dxi0, dxj0, dxk0);

    double dxi02 = dxi0*dxi0;
    double dxj02 = dxj0*dxj0;
    double dxk02 = dxk0*dxk0;

    closestBoundary(y, 1)->closestImage(y, x, z, yp, xp, zp, dxj1, dxi1, dxk1);

    double dx2 = dxi1*dxi1;
    double dy2 = dxj1*dxj1;
    double dz2 = dxk1*dxk1;

    if (dxi02 < dx2)
    {
        dx2 = dxi02;
    }

    if (dxj02 < dy2)
    {
        dy2 = dxj02;
    }

    if (dxk02 < dz2)
    {
        dz2 = dxk02;
    }


    double dxk2 = z - zp;
    double dxk22 = dxk2*dxk2;

    if (dxk22 < dz2)
    {
        dz2 = dxk22;
    }

    return dx2 + dy2 + dz2;


}

double SOSSolver::absSquareDistance(const uint x, const uint y, const int z,
                                    const double xp, const double yp, const double zp) const
{
    double dx = x - xp;
    double dy = y - yp;
    double dz = z - zp;

    return dx*dx + dy*dy + dz*dz;
}

void SOSSolver::execute()
{
    KMCSolver::execute();

    m_changedSurfaceSites.clear();
}


void SOSSolver::initialize()
{
    if (!m_heights_set)
    {

        m_heights.zeros();

        int left;
        int bottom;

        m_averageHeight = 0;

        for (uint x = 0; x < m_length; ++x)
        {
            for (uint y = 0; y < m_width; ++y)
            {
                left = leftSite(x, y, m_heights(x, y));
                bottom = bottomSite(x, y, m_heights(x, y));

                if (isOutsideBoxSingle(left, 0) || isOutsideBoxSingle(bottom, 1))
                {
                    continue;
                }

                m_heights(x, y) = (m_heights(left, y) + m_heights(x, bottom))/2 + round(-1 + 2*rng.uniform());
            }
        }
    }

    m_averageHeight = arma::accu(m_heights)/double(area());

    for (uint x = 0; x < m_length; ++x)
    {
        for (uint y = 0; y < m_width; ++y)
        {
            setNNeighbors(x, y);
        }
    }

    for (HeightConnecter *connecter : m_heightConnecters)
    {
        connecter->setupInitialConditions();
    }

    KMCSolver::initialize();

}
