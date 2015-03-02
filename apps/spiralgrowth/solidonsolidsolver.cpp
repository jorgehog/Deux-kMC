#include "solidonsolidsolver.h"
#include "solidonsolidreaction.h"
#include "pressurewall.h"


SolidOnSolidSolver::SolidOnSolidSolver(const uint length,
                                       const uint width,
                                       const double alpha,
                                       const double mu,
                                       const bool shadowing) :
    KMCSolver(),
    m_dim((( length == 1 ) || ( width == 1 ) ) ? 1 : 2),
    m_length(length),
    m_width(width),
    m_alpha(alpha),
    m_mu(mu),
    m_shadowing(shadowing),
    m_heights(length, width, fill::zeros),
    m_nNeighbors(length, width),
    m_siteReactions(length, width)
{
    for (uint x = 0; x < length; ++x)
    {
        for (uint y = 0; y < width; ++y)
        {
            m_heights(x, y) = (m_heights(leftSite(x), y) + m_heights(x, bottomSite(y)))/2 + round(-1 + 2*rng.uniform());
            m_siteReactions(x, y) = new DiffusionDeposition(x, y, *this);
        }
    }

    for (uint x = 0; x < length; ++x)
    {
        for (uint y = 0; y < width; ++y)
        {
            setNNeighbors(x, y);
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

    setNNeighbors(x, y);
    setNNeighbors(leftSite(x), y);
    setNNeighbors(rightSite(x), y);
    setNNeighbors(x, bottomSite(y));
    setNNeighbors(x, topSite(y));
}

void SolidOnSolidSolver::setNNeighbors(const uint x, const uint y)
{
    m_nNeighbors(x, y) = calculateNNeighbors(x, y);
}


void SolidOnSolidSolver::setPressureWallEvent(PressureWall &pressureWallEvent)
{
    m_pressureWallEvent = &pressureWallEvent;
}

double SolidOnSolidSolver::localPressure(const uint x, const uint y) const
{
    if (!m_pressureWallEvent->hasStarted())
    {
        return 0;
    }

    return m_pressureWallEvent->localPressure(x, y);
}

uint SolidOnSolidSolver::calculateNNeighbors(const uint x, const uint y) const
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
    return 3 - n/2;
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
                _reaction.setDiffusionRate(_reaction.diffusionRate()*expFac);
                _reaction.changeRate(_reaction.diffusionRate() + _reaction.depositionRate());
            }
        }
    }

    m_mu = mu;

}

uint SolidOnSolidSolver::numberOfReactions() const
{
    return m_siteReactions.size();
}

Reaction *SolidOnSolidSolver::getReaction(const uint n) const
{
    return m_siteReactions(n);
}
