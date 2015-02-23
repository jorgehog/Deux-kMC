#pragma once

#include <kMC>
#include <armadillo>

using namespace kMC;
using namespace arma;

class DiffusionDeposition;
class PressureWall;

class SolidOnSolidSolver : public KMCSolver
{
public:
    SolidOnSolidSolver(PressureWall &pressureWallEvent,
                       const uint length,
                       const uint width,
                       const double alpha,
                       const double mu,
                       const bool shadowing);

    ~SolidOnSolidSolver();

    void registerHeightChange(const uint x, const uint y, const int value)
    {
        m_heights(x, y) += value;
    }

    const uint &length() const
    {
        return m_length;
    }

    const uint &width() const
    {
        return m_width;
    }

    uint area() const
    {
        return m_length*m_width;
    }

    const uint &dim() const
    {
        return m_dim;
    }

    const double &alpha() const
    {
        return m_alpha;
    }

    const double &mu() const
    {
        return m_mu;
    }

    const bool &shadowing() const
    {
        return m_shadowing;
    }

    const imat &heights() const
    {
        return m_heights;
    }

    const int &height(const uint x, const uint y) const
    {
        return m_heights(x, y);
    }

    PressureWall &pressureWallEvent() const
    {
        return m_pressureWallEvent;
    }

    double localPressure(const uint x, const uint y) const;

    uint nNeighbors(const uint x, const uint y) const;

    uint topSite(const uint site, const uint n = 1) const;

    uint bottomSite(const uint site, const uint n = 1) const;

    uint leftSite(const uint site, const uint n = 1) const;

    uint rightSite(const uint site, const uint n = 1) const;

    DiffusionDeposition &reaction(const uint x, const uint y) const
    {
        return *m_siteReactions(x, y);
    }

    double shadowScale(const double n) const;

    void setMu(const double mu)
    {
        m_mu = mu;
    }

private:

    const uint m_dim;

    PressureWall &m_pressureWallEvent;

    const uint m_length;
    const uint m_width;

    const double m_alpha;
    double m_mu;

    const bool m_shadowing;

    imat m_heights;

    field<DiffusionDeposition*> m_siteReactions;

    // KMCSolver interface
public:
    uint numberOfReactions() const;
    Reaction *getReaction(const uint n) const;
};
