#pragma once

#include "../kmcsolver/kmcsolver.h"
#include <armadillo>

using namespace kMC;
using namespace arma;

class DissolutionDeposition;
class ConfiningSurface;
class Diffusion;

class SOSSolver : public KMCSolver
{
public:
    SOSSolver(const uint length,
              const uint width,
              const double alpha,
              const double gamma,
              const Boundary* xBoundary,
              const Boundary* yBoundary);

    virtual ~SOSSolver();

    void registerHeightChange(const uint x, const uint y, const int value);

    void setNNeighbors(const uint x, const uint y);

    void setHeight(const uint x, const uint y, const int value);

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

    uint dim() const
    {
        return surfaceDim() + 1;
    }

    const uint &surfaceDim() const
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

    double concentration() const
    {
        return exp(gamma() - 2*alpha());
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

    void setDiffusionEvent(Diffusion &diffusionEvent);

    ConfiningSurface &confiningSurfaceEvent() const
    {
        return *m_confiningSurfaceEvent;
    }

    Diffusion &diffusionEvent() const
    {
        return *m_diffusionEvent;
    }

    double volume() const;

    double confinementEnergy(const uint x, const uint y) const;

    double depositionRate(const uint x, const uint y) const;

    uint calculateNNeighbors(const uint x, const uint y) const;

    uint nSurroundingSolutionSites(const uint x, const uint y) const;

    void getSolutionSite(const uint x, const uint y,
                         int &dx, int &dy, int &dz,
                         const uint siteNumber) const;

    int topSite(const uint site, const uint n = 1) const;

    int bottomSite(const uint site, const uint n = 1) const;

    int leftSite(const uint site, const uint n = 1) const;

    int rightSite(const uint site, const uint n = 1) const;

    void findConnections(const uint x,
                         const uint y,
                         bool &connectedLeft,
                         bool &connectedRight,
                         bool &connectedBottom,
                         bool &connectedTop) const;

    uint span() const;

    bool isBlockedPosition(const double x, const double y, const double z) const;

    DissolutionDeposition &surfaceReaction(const uint x, const uint y) const
    {
        return *m_siteReactions(x, y);
    }

    void setMu(const double gamma);

private:

    const uint m_dim;

    ConfiningSurface *m_confiningSurfaceEvent;
    Diffusion *m_diffusionEvent;

    const uint m_length;
    const uint m_width;

    const double m_alpha;
    double m_mu;

    imat m_heights;
    umat m_nNeighbors;

    field<DissolutionDeposition*> m_siteReactions;

    std::vector<DissolutionDeposition*> m_affectedReactions;

    // KMCSolver interface
private:
    void initializeSolver();

};
