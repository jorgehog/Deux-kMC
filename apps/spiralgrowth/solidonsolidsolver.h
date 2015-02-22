#pragma once

#include <kMC>
#include <armadillo>

using namespace kMC;
using namespace arma;

class SolidOnSolidSolver : public KMCSolver
{
public:
    SolidOnSolidSolver(const uint length, const uint width, const double alpha);
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

    const imat &heights() const
    {
        return m_heights;
    }

    const int &height(const uint x, const uint y) const
    {
        return m_heights(x, y);
    }

    uint nNeighbors(const uint x, const uint y) const;

    uint topSite(const uint site, const uint n = 1) const;

    uint bottomSite(const uint site, const uint n = 1) const;

    uint leftSite(const uint site, const uint n = 1) const;

    uint rightSite(const uint site, const uint n = 1) const;

    Reaction &reaction(const uint x, const uint y) const
    {
        return *m_siteReactions(x, y);
    }

private:

    const uint m_dim;

    const uint m_length;
    const uint m_width;

    const double m_alpha;

    imat m_heights;

    field<Reaction*> m_siteReactions;

    // KMCSolver interface
public:
    uint numberOfReactions() const;

private:
    Reaction *getReaction(const uint n) const;
};
