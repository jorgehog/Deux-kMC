#pragma once

#include "../kmcsolver/kmcsolver.h"
#include <armadillo>

using namespace kMC;
using namespace arma;

class DiffusionDeposition;
class ConfiningSurface;
class CavityDiffusion;

class SolidOnSolidSolver : public KMCSolver
{
public:
    SolidOnSolidSolver(const uint length,
                       const uint width,
                       const Boundary* xBoundary,
                       const Boundary* yBoundary,
                       const double alpha,
                       const double gamma);

    ~SolidOnSolidSolver();

    void registerHeightChange(const uint x, const uint y, const int value);

    void setNNeighbors(const uint x, const uint y);

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

    const double &gamma() const
    {
        return m_mu;
    }

    const imat &heights() const
    {
        return m_heights;
    }

    const int &height(const uint x, const uint y) const
    {
        return m_heights(x, y);
    }

    const uint &nNeighbors(const uint x, const uint y) const
    {
        BADAssEqual(m_nNeighbors(x, y), calculateNNeighbors(x, y));

        return m_nNeighbors(x, y);
    }

    void setConfiningSurfaceEvent(ConfiningSurface &confiningSurfaceEvent);

    void setDiffusionEvent(CavityDiffusion &diffusionEvent);

    ConfiningSurface &confiningSurfaceEvent() const
    {
        return *m_confiningSurfaceEvent;
    }

    CavityDiffusion &diffusionEvent() const
    {
        return *m_diffusionEvent;
    }

    double confinementEnergy(const uint x, const uint y) const;

    double localSurfaceSupersaturation(const uint x, const uint y) const;

    uint calculateNNeighbors(const uint x, const uint y) const;

    int topSite(const uint site, const uint n = 1) const;

    int bottomSite(const uint site, const uint n = 1) const;

    int leftSite(const uint site, const uint n = 1) const;

    int rightSite(const uint site, const uint n = 1) const;

    uint span() const;

    DiffusionDeposition &reaction(const uint x, const uint y) const
    {
        return *m_siteReactions(x, y);
    }

    void setMu(const double gamma);

private:

    const uint m_dim;

    ConfiningSurface *m_confiningSurfaceEvent;
    CavityDiffusion *m_diffusionEvent;

    const uint m_length;
    const uint m_width;

    const double m_alpha;
    double m_mu;

    imat m_heights;
    umat m_nNeighbors;

    field<DiffusionDeposition*> m_siteReactions;

    void initializeSolver();

    std::vector<DiffusionDeposition*> m_affectedReactions;


    // KMCSolver interface
public:
    uint numberOfReactions() const;
    Reaction *getReaction(const uint n) const;
};
