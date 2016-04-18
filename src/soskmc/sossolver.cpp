#include "sossolver.h"
#include "dissolutiondeposition.h"
#include "concentrationboundaryreaction.h"
#include "Events/confiningsurface/confiningsurface.h"
#include "Events/diffusion/constantconcentration.h"
#include "../kmcsolver/boundary/boundary.h"


SOSSolver::SOSSolver(const uint length,
                     const uint width,
                     const double alpha,
                     const double gamma,
                     const bool surfaceDiffusion) :
    KMCSolver(),
    Subject(),
    m_surfaceDiffusion(surfaceDiffusion),
    m_heightsSet(false),
    m_surfaceDim((( length == 1 ) || ( width == 1 ) ) ? 1 : 2),
    m_confiningSurfaceEvent(nullptr),
    m_diffusionEvent(nullptr),
    m_length(length),
    m_width(width),
    m_alpha(alpha),
    m_c0(exp(-(int)dim()*alpha)),
    m_gamma(gamma),
    m_concentration(exp(gamma-(int)dim()*alpha)),
    m_expGamma(exp(gamma)),
    m_concentrationIsZero(false),
    m_heights(length, width, fill::zeros),
    m_nNeighbors(length, width),
    m_surfaceReactions(length, width)
{
    for (uint x = 0; x < length; ++x)
    {
        for (uint y = 0; y < width; ++y)
        {
            m_surfaceReactions(x, y) = new SurfaceReaction(x, y, *this);
            addReaction(m_surfaceReactions(x, y));
        }
    }
}

SOSSolver::~SOSSolver()
{
    for (uint x = 0; x < length(); ++x)
    {
        for (uint y = 0; y < width(); ++y)
        {
            delete m_surfaceReactions(x, y);
        }
    }

    for (ConcentrationBoundaryReaction *r : m_concentrationBoundaryReactions)
    {
        delete r;
    }

    m_surfaceReactions.clear();
    m_concentrationBoundaryReactions.clear();
}

void SOSSolver::registerHeightChange(const uint x, const uint y, const int value)
{
    BADAssBool(!isOutsideBox(x, y));
    BADAssEqual(abs(value), 1);

    m_heights(x, y) += value;

    registerChangedSite(x, y);

    m_averageHeight += value/double(area());
    BADAssClose(averageHeight(), accu(heights())/double(area()), 1E-3);

    registerChangedAround(x, y);

    m_currentSurfaceChange.x = x;
    m_currentSurfaceChange.y = y;
    m_currentSurfaceChange.value = value;
    m_currentSurfaceChange.type = ChangeTypes::Single;

    notifyObservers(Subjects::SOLVER);
}

void SOSSolver::registerSurfaceTransition(const uint x0, const uint y0, const int x1, const int y1)
{
    BADAssBool(!isOutsideBox(x0, y0), "invalid site.");
    BADAssBool(isSurfaceSite(x1, y1, height(x0, y0)), "ERRAH.");

    m_heights(x0, y0) -= 1;
    registerChangedSite(x0, y0);
    registerChangedAround(x0, y0);

    if (!isOutsideBox(x1, y1))
    {
        m_heights(x1, y1) += 1;
        registerChangedSite(x1, y1);
        registerChangedAround(x1, y1);
    }

    else
    {
        m_averageHeight -= 1./area();
    }

    m_currentSurfaceChange.x = x0;
    m_currentSurfaceChange.y = y0;
    m_currentSurfaceChange.x1 = x1;
    m_currentSurfaceChange.y1 = y1;
    m_currentSurfaceChange.type = ChangeTypes::Double;

    notifyObservers(Subjects::SOLVER);
}

void SOSSolver::registerChangedSite(const uint x, const uint y)
{
    m_changedSurfaceSites.insert(make_pair(x, y));
    setNNeighbors(x, y);
    registerAffectedReaction(&surfaceReaction(x, y));
}

void SOSSolver::registerChangedAround(const uint x, const uint y)
{
    const uint left = leftSite(x, y, m_heights(x, y));
    const uint right = rightSite(x, y, m_heights(x, y));
    const uint bottom = bottomSite(x, y, m_heights(x, y));
    const uint top = topSite(x, y, m_heights(x, y));

    if (!isOutsideBoxSingle(left, 0))
    {
        setNNeighbors(left, y);
        registerAffectedReaction(&surfaceReaction(left, y));
    }

    if (!isOutsideBoxSingle(right, 0))
    {
        setNNeighbors(right, y);
        registerAffectedReaction(&surfaceReaction(right, y));
    }

    if (!isOutsideBoxSingle(top, 1))
    {
        setNNeighbors(x, top);
        registerAffectedReaction(&surfaceReaction(x, top));
    }

    if (!isOutsideBoxSingle(bottom, 1))
    {
        setNNeighbors(x, bottom);
        registerAffectedReaction(&surfaceReaction(x, bottom));
    }
}

void SOSSolver::setNNeighbors(const uint x, const uint y)
{
    m_nNeighbors(x, y) = calculateNNeighbors(x, y);
}

void SOSSolver::setHeight(const uint x, const uint y, const int value, const bool iteratively)
{
    if (!m_heightsSet)
    {
        m_heightsSet = true;
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

void SOSSolver::setHeights(const imat &newheights, const bool iteratively)
{
    for (uint x = 0; x < length(); ++x)
    {
        for (uint y = 0; y < width(); ++y)
        {
            setHeight(x, y, newheights(x, y), iteratively);
        }
    }
}

bool SOSSolver::depositionIsAvailable(const uint x, const uint y) const
{
    if (m_confiningSurfaceEvent->hasSurface())
    {
        if (m_heights(x, y) + 1 > m_confiningSurfaceEvent->height() - 1)
        {
            return false;
        }
    }

    return true;
}

void SOSSolver::setConfiningSurfaceEvent(ConfiningSurface &confiningSurfaceEvent)
{
    m_confiningSurfaceEvent = &confiningSurfaceEvent;
    confiningSurfaceEvent.registerObserver(this);

    if (diffusionEventIsSet())
    {
        confiningSurfaceEvent.registerObserver(m_diffusionEvent);
        registerObserver(m_diffusionEvent);
        registerObserver(m_confiningSurfaceEvent);
    }
}

void SOSSolver::setDiffusionEvent(Diffusion &diffusionEvent)
{
    m_diffusionEvent = &diffusionEvent;

    if (confiningSurfaceIsSet())
    {
        confiningSurfaceEvent().registerObserver(&diffusionEvent);
        registerObserver(&diffusionEvent);
        registerObserver(m_confiningSurfaceEvent);
    }
}

double SOSSolver::volume() const
{
    //sum (h_l - h_i - 1)
    return (m_confiningSurfaceEvent->height() - 1)*area() - arma::accu(m_heights);
}

void SOSSolver::calculateHeightDependentValues()
{
    m_averageHeight = arma::accu(m_heights)/double(area());

    for (uint x = 0; x < m_length; ++x)
    {
        for (uint y = 0; y < m_width; ++y)
        {
            setNNeighbors(x, y);
        }
    }

    //initialize surface on top of the current surface as default.
    if (confiningSurfaceEvent().hasSurface())
    {
        const int maxHeight = m_heights.max();
        if (confiningSurfaceEvent().height() == 0 && maxHeight > 0)
        {
            confiningSurfaceEvent().setHeight(maxHeight + 2);
        }
    }

    initializeObservers(Subjects::SOLVER);
    confiningSurfaceEvent().initializeObservers(Subjects::CONFININGSURFACE);
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

uint SOSSolver::numberOfSurroundingSites(const uint x, const uint y)
{
    bool connectedLeft, connectedRight,
            connectedTop, connectedBottom;

    findConnections(x, y,
                    connectedLeft,
                    connectedRight,
                    connectedBottom,
                    connectedTop,
                    false);

    const int &h = height(x, y);

    bool connectedAbove = diffusionEvent().isBlockedPosition(x, y, h+1);

    if (confiningSurfaceEvent().hasSurface())
    {
        connectedAbove = connectedAbove || isBlockedPosition(x, y, h+1);
    }

    uint n = 5;

    if (connectedLeft)
    {
        n--;
    }

    if (connectedRight)
    {
        n--;
    }

    if (connectedBottom)
    {
        n--;
    }

    if (connectedTop)
    {
        n--;
    }

    if (connectedAbove)
    {
        n--;
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

void SOSSolver::getRandomSolutionSite(const uint x, const uint y, const int height, int &dx, int &dy, int &dz) const
{
    const uint nSites = numberOfSurroundingSolutionSites(x, y, height);
    const uint randomSite = rng.uniform()*nSites;

    BADAss(nSites, !=, 0u);

    getSolutionSite(x, y, height, dx, dy, dz, randomSite);

}

int SOSSolver::topSite(const uint x, const uint y, const int z, const uint n) const
{
    return boundary(1, 1)->transformLatticeCoordinate(y + n, x, z);
}

int SOSSolver::bottomSite(const uint x, const uint y, const int z, const uint n) const
{
    return boundary(1, 0)->transformLatticeCoordinate((int)y - (int)n, x, z);
}

int SOSSolver::leftSite(const uint x, const uint y, const int z, const uint n) const
{
    return boundary(0, 0)->transformLatticeCoordinate((int)x - (int)n, y, z);
}

int SOSSolver::rightSite(const uint x, const uint y, const int z, const uint n) const
{
    return boundary(0, 1)->transformLatticeCoordinate(x + n, y, z);
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

    if (!boundary(dim, orientation)->isBlockedLattice(xNeighbor, y, h))
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

const Boundary *SOSSolver::getBoundaryFromLoc(const int x, const int y) const
{
    if (x < 0)
    {
        return boundary(0, 0);
    }

    else if (y < 0)
    {
        return boundary(1, 0);
    }

    else if (uint(x) >= length())
    {
        return boundary(0, 1);
    }

    else
    {
        return boundary(1, 1);
    }
}

void SOSSolver::boundaryLatticeTransform(int &xTrans, int &yTrans, const int x, const int y, const int z) const
{
    uint xOrientation;
    uint yOrientation;

    if (x < signed(length()/2))
    {
        xOrientation = 0;
    }

    else
    {
        xOrientation = 1;
    }

    if (y < signed(width()/2))
    {
        yOrientation = 0;
    }

    else
    {
        yOrientation = 1;
    }

    xTrans = boundary(0, xOrientation)->transformLatticeCoordinate(x, y, z);
    yTrans = boundary(1, yOrientation)->transformLatticeCoordinate(y, x, z);
}

void SOSSolver::boundaryContinousTransform(double &xTrans, double &yTrans, const double x, const double y, const double z) const
{
    uint xOrientation;
    uint yOrientation;

    if (x < length()/2)
    {
        xOrientation = 0;
    }

    else
    {
        xOrientation = 1;
    }

    if (y < width()/2)
    {
        yOrientation = 0;
    }

    else
    {
        yOrientation = 1;
    }

    xTrans = boundary(0, xOrientation)->transformContinousCoordinate(x, y, z);
    yTrans = boundary(1, yOrientation)->transformContinousCoordinate(y, x, z);
}

int SOSSolver::boundaryLatticeTransformSingle(const int x, const int y, const int z, uint dim, const int shift) const
{
    uint orientation;

    if (dim == 0)
    {
        if (x < signed(length()/2))
        {
            orientation = 0;
        }

        else
        {
            orientation = 1;
        }

        return boundary(0, orientation)->transformLatticeCoordinate(x + shift, y, z);
    }

    else
    {
        if (y < signed(width()/2))
        {
            orientation = 0;
        }

        else
        {
            orientation = 1;
        }

        return boundary(0, orientation)->transformLatticeCoordinate(y + shift, x, z);
    }

}

double SOSSolver::boundaryContinousTransformSingle(const double x, const double y, const double z, uint dim, const double shift) const
{
    uint orientation;

    if (dim == 2)
    {
        return z + shift;
    }

    else if (dim == 0)
    {
        if (x < length()/2)
        {
            orientation = 0;
        }

        else
        {
            orientation = 1;
        }

        return boundary(0, orientation)->transformContinousCoordinate(x + shift, y, z);
    }

    else
    {
        if (y < width()/2)
        {
            orientation = 0;
        }

        else
        {
            orientation = 1;
        }

        return boundary(0, orientation)->transformContinousCoordinate(y + shift, x, z);
    }
}

//uint SOSSolver::boundaryOrientation(const double x, const uint dim) const
//{
//    const uint lx = dim == 0 ? length() : width();

//    return x >= lx/2 ? 1 : 0;
//}

//void SOSSolver::boundaryTransformXY(double &xTrans, double &yTrans, const double x, const double y, const double z) const
//{
//    uint xOrientation;
//    uint yOrientation;

//    if (x < length()/2)
//    {
//        xOrientation = 0;
//    }

//    else
//    {
//        xOrientation = 1;
//    }

//    if (y < width()/2)
//    {
//        yOrientation = 0;
//    }

//    else
//    {
//        yOrientation = 1;
//    }

//    xTrans = boundary(0, xOrientation)->transformCoordinate(x, y, z);
//    yTrans = boundary(1, yOrientation)->transformCoordinate(x, y, z);

//}

//double SOSSolver::boundaryTransform(const double x, const double y, const double z, const uint dim) const
//{
//    uint orientation;

//    if (dim == 0)
//    {
//        if (x < length()/2)
//        {
//            orientation = 0;
//        }

//        else
//        {
//            orientation = 1;
//        }

//        return boundary(0, orientation)->transformCoordinate(x, y, z);
//    }

//    else if (dim == 1)
//    {
//        if (y < length()/2)
//        {
//            orientation = 0;
//        }

//        else
//        {
//            orientation = 1;
//        }

//        return boundary(1, orientation)->transformCoordinate(y, x, z);
//    }

//    else
//    {
//        return z;
//    }

//}

//double SOSSolver::boundaryTransform(const double x, const double y, const double z, const double dxi, const uint dim) const
//{

//    if (dim == 0)
//    {
//        return closestBoundary(x, dim)->transformCoordinate(x + dxi, y, z);
//    }

//    else if (dim == 1)
//    {
//        return closestBoundary(x, dim)->transformCoordinate(y + dxi, x, z);
//    }

//    else
//    {
//        return z + dxi;
//    }

//}

//const Boundary *SOSSolver::closestBoundary(const double x, const uint dim) const
//{
//    return boundary(dim, boundaryOrientation(x, dim));
//}

void SOSSolver::addConcentrationBoundary(const uint dim,
                                         const Boundary::orientations orientation,
                                         const double omega)
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

    ConcentrationBoundaryReaction *concReaction = new ConcentrationBoundaryReaction(dim, orientationInt, *this, omega);

    m_concentrationBoundaryReactions.push_back(concReaction);

    addReaction(concReaction);

    registerObserver(concReaction);
    confiningSurfaceEvent().registerObserver(concReaction);
}

bool SOSSolver::isBlockedPosition(const double x, const double y, const double z) const
{
    //center of conf surf is at confSE().height(), such that it extends to confSE().height() - 0.5
    //which makes a particle of radius 0.5 collide if z is larger than cse.h() - 1
    bool isOutSideBox_z;
    if (confiningSurfaceEvent().hasSurface())
    {
        isOutSideBox_z = (z > confiningSurfaceEvent().height() - 1);
    }
    else
    {
        isOutSideBox_z = false;
    }

    if (isOutSideBox_z)
    {
        return true;
    }

    const uint X = uint(round(x));
    const uint Y = uint(round(y));

    BADAssBool(!isOutsideBox(X, Y));

    //center of surface particle is at h=height(X, Y) and surface particles extend to h+0.5
    //such that particles of radius 0.5 collide when z is lower than h - 1
    return z < height(X, Y) + 1;
}

bool SOSSolver::isOutsideBoxContinuous(const double x, const double y) const
{
    //issue here when x = 19.7 it is outside, but not really. it is at 0.2.
    //if we let 19.7 pass it will round to 20. Need periodic to
    //transform 19.7 into -0.2 so it rounds up to 0

    //to enable round arond 0 and length to have same room as othes we let particles
    //go from -0.5 to l-0.5. Boundaries should not suggest these moves
    //if they are illegal.
    const bool isOutSideBox_x = (x < -0.5) || (x > length() - 0.5);

    const bool isOutSideBox_y = (y < -0.5) || (y > width() - 0.5);

    return (isOutSideBox_x || isOutSideBox_y);
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
    if (x < 0 || y < 0)
    {
        return true;
    }

    if (uint(x) >= length() || uint(y) >= width())
    {
        return true;
    }

    return false;
}

bool SOSSolver::isSurfaceSite(const int x, const int y, const int z) const
{
    if (isOutsideBox(x, y))
    {
        const Boundary *b = getBoundaryFromLoc(x, y);

        return b->isBlockedLattice(x, y, z-1) && (!b->isBlockedLattice(x, y, z));
    }

    return height(x, y) == z - 1;
}

void SOSSolver::setGamma(const double gamma)
{
    const double gammaPrev = m_gamma;

    m_gamma = gamma;
    m_concentration = exp(m_gamma - (int)dim()*alpha());
    m_expGamma = exp(m_gamma);

    if (hasStarted())
    {

        if (concentrationIsZero())
        {
            m_concentrationIsZero = false;
            recalculateAllRates();
        }

        else
        {

            double expFac = exp(gammaPrev - gamma);
            BADAss(expFac, ==, expFac);

            for (uint x = 0; x < length(); ++x)
            {
                for (uint y = 0; y < width(); ++y)
                {
                    SurfaceReaction &_reaction = surfaceReaction(x, y);
                    _reaction.setEscapeRate(_reaction.escapeRate()*expFac);
                }
            }
        }
    }

}

void SOSSolver::setConcentration(const double concentration)
{
    double gamma = log(concentration) + dim()*alpha();

    setGamma(gamma);
}

void SOSSolver::shiftConcentration(const double dc)
{
    setConcentration(concentration() + dc);
}

double SOSSolver::closestSquareDistance(const uint x, const uint y, const int z,
                                        const double xp, const double yp, const double zp) const
{
    uint xOrientation;
    uint yOrientation;

    if (x < length()/2)
    {
        xOrientation = 0;
    }

    else
    {
        xOrientation = 1;
    }

    if (y < width()/2)
    {
        yOrientation = 0;
    }

    else
    {
        yOrientation = 1;
    }

    double dxi0, dxj0, dxk0;
    double dxi1, dxj1, dxk1;

    boundary(0, xOrientation)->closestImage(x, y, z, xp, yp, zp, dxi0, dxj0, dxk0);

    double dxi02 = dxi0*dxi0;
    double dxj02 = dxj0*dxj0;
    double dxk02 = dxk0*dxk0;

    boundary(1, yOrientation)->closestImage(y, x, z, yp, xp, zp, dxj1, dxi1, dxk1);

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

void SOSSolver::setZeroConcentration()
{
    m_concentrationIsZero = true;

    m_gamma = std::numeric_limits<double>::min();
    m_concentration = 0;
    m_expGamma = 0;

    if (hasStarted())
    {
        recalculateAllRates();
    }
}

void SOSSolver::freezeSurfaceParticle(const uint x, const uint y)
{
    m_surfaceReactions(x, y)->freeze();
}

void SOSSolver::execute()
{
    KMCSolver::execute();

    m_changedSurfaceSites.clear();
}


void SOSSolver::initialize()
{
    if (!m_heightsSet)
    {
        m_heights.zeros();
    }

    calculateHeightDependentValues();

    KMCSolver::initialize();

}

double SOSSolver::timeUnit() const
{
    //original rate expression is exp(-alpha*n)/c0

    //divide by nothing
    if (concentrationIsZero())
    {
        return 1.;
    }

    //divide rates by c gives multiplication of timestep by same factor
    else
    {
        return 1./concentration();
    }
}

void SOSSolver::initializeObserver(const Subjects &subject)
{
    (void) subject;
}

void SOSSolver::notifyObserver(const Subjects &subject)
{
    (void) subject;

    const double &h = confiningSurfaceEvent().height();
    const int surfaceContactHeight = floor(h) - 1;

    for (uint x = 0; x < length(); ++x)
    {
        for (uint y = 0; y < width(); ++y)
        {
            if (height(x, y) == surfaceContactHeight || height(x, y) == surfaceContactHeight - 1)
            {
                registerAffectedReaction(m_surfaceReactions(x, y));
            }
        }
    }

}
