#include "solidonsolidsolver.h"
#include "solidonsolidreaction.h"
#include "pressurewall.h"


SolidOnSolidSolver::SolidOnSolidSolver(PressureWall &pressureWallEvent, const uint length,
                                       const uint width,
                                       const double alpha,
                                       const double mu,
                                       const bool shadowing) :
    KMCSolver(),
    m_dim((( length == 1 ) || ( width == 1 ) ) ? 1 : 2),
    m_pressureWallEvent(pressureWallEvent),
    m_length(length),
    m_width(width),
    m_alpha(alpha),
    m_mu(mu),
    m_shadowing(shadowing),
    m_heights(length, width),
    m_siteReactions(length, width)
{
    m_pressureWallEvent.setDependency(*this);

    for (uint x = 0; x < length; ++x)
    {
        for (uint y = 0; y < width; ++y)
        {
            m_heights(x, y) = (m_heights(leftSite(x), y) + m_heights(x, bottomSite(y)))/2 + round(-1 + 2*rng.uniform());
            m_siteReactions(x, y) = new DiffusionDeposition(x, y, *this);
        }
    }

    m_pressureWallEvent.setupInitialConditions();

}

SolidOnSolidSolver::~SolidOnSolidSolver()
{
    for (uint x = 0; x < length(); ++x)
    {
        for (uint y = 0; y < width(); ++y)
        {
            delete m_siteReactions(x);
        }
    }

    m_siteReactions.clear();
}

double SolidOnSolidSolver::localPressure(const uint x, const uint y) const
{
    if (!m_pressureWallEvent.hasStarted())
    {
        return 0;
    }

    return m_pressureWallEvent.localPressure(x, y);
}

uint SolidOnSolidSolver::nNeighbors(const uint x, const uint y) const
{

    uint n = 1;

    const int &h = m_heights(x, y);
    const int &hLeft = m_heights(leftSite(x), y);
    const int &hRight = m_heights(rightSite(x), y);
    const int &hTop = m_heights(x, topSite(y));
    const int &hBottom = m_heights(x, bottomSite(y));


    bool connectedLeft = hLeft >= h;
    bool connectedRight = hRight >= h;
    bool connectedTop = hTop >= h;
    bool connectedBottom = hBottom >= h;

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

uint SolidOnSolidSolver::topSite(const uint site, const uint n) const
{
    return (site + n)%width();
}

uint SolidOnSolidSolver::bottomSite(const uint site, const uint n) const
{
    return (site + width() - n)%width();
}

uint SolidOnSolidSolver::leftSite(const uint site, const uint n) const
{
    return (site + length() - n)%length();
}

uint SolidOnSolidSolver::rightSite(const uint site, const uint n) const
{
    return (site + n)%length();
}

double SolidOnSolidSolver::shadowScale(const double n) const
{
    return 2*(1 + m_dim) - n;
}

uint SolidOnSolidSolver::numberOfReactions() const
{
    return m_siteReactions.size();
}

Reaction *SolidOnSolidSolver::getReaction(const uint n) const
{
    return m_siteReactions(n);
}
