#pragma once

#include "../kmcsolver/kmcsolver.h"
#include "../kmcsolver/boundary/boundary.h"
#include <armadillo>

#include <unordered_set>
#include <boost/functional/hash.hpp>

using namespace kMC;
using namespace arma;
using std::pair;
using std::unordered_set;

class DissolutionDeposition;
class ConcentrationBoundaryReaction;
class ConfiningSurface;
class Diffusion;

class SOSSolver : public KMCSolver
{
    typedef pair<uint, uint> pair_type;
    typedef unordered_set<pair_type, boost::hash<pair_type>> set_type;

public:
    SOSSolver(const uint length,
              const uint width,
              const double alpha,
              const double gamma,
              vector<vector<const Boundary *> > boundaries);

    virtual ~SOSSolver();

    void registerHeightChange(const uint x, const uint y, const int value);

    void setNNeighbors(const uint x, const uint y);

    void setHeight(const uint x, const uint y, const int value, const bool iteratively = true);

    void setHeights(const imat &heights, const bool iteratively = true);

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

    uint calculateNNeighbors(const uint x, const uint y, const int h) const;

    uint calculateNNeighbors(const uint x, const uint y) const
    {
        return calculateNNeighbors(x, y, height(x, y));
    }

    uint numberOfSurroundingSolutionSites(const uint x, const uint y, const int h) const;

    uint numberOfSurroundingSolutionSites(const uint x, const uint y) const
    {
        return numberOfSurroundingSolutionSites(x, y, height(x, y));
    }

    void getSolutionSite(const uint x, const uint y, const int height,
                         int &dx, int &dy, int &dz,
                         const uint siteNumber) const;

    void getSolutionSite(const uint x, const uint y,
                         int &dx, int &dy, int &dz,
                         const uint siteNumber) const
    {
        getSolutionSite(x, y, height(x, y), dx, dy, dz, siteNumber);
    }

    int topSite(const uint site, const uint n = 1) const;

    int bottomSite(const uint site, const uint n = 1) const;

    int leftSite(const uint site, const uint n = 1) const;

    int rightSite(const uint site, const uint n = 1) const;

    void findConnections(const uint x,
                         const uint y,
                         const int h,
                         bool &connectedLeft,
                         bool &connectedRight,
                         bool &connectedBottom,
                         bool &connectedTop,
                         bool onlySurface = true) const;

    void findConnections(const uint x,
                         const uint y,
                         bool &connectedLeft,
                         bool &connectedRight,
                         bool &connectedBottom,
                         bool &connectedTop,
                         bool onlySurface = true) const
    {
        return findConnections(x, y, height(x, y),
                               connectedLeft,
                               connectedRight,
                               connectedBottom,
                               connectedTop,
                               onlySurface);
    }

    bool findSingleConnection(const int neighbor, const uint dim, const uint orientation, const uint x, const int h, bool onlySurface) const;

    uint span() const;

    uint xBoundaryOrientation(const uint x) const;

    uint yBoundaryOrientation(const uint y) const;

    uint boundaryOrientation(const uint x, const uint lx) const;

    int boundaryTransform(const uint x, const int dx, const uint dim) const;

    void addConcentrationBoundary(const uint dim, const Boundary::orientations orientation);

    bool isBlockedPosition(const double x, const double y, const double z) const;

    bool isOutsideBoxSingle(const int x, const uint dim) const;

    bool isOutsideBox(const int x, const int y) const;

    bool isSurfaceSite(const uint x, const uint y, const int z) const;

    DissolutionDeposition &surfaceReaction(const uint x, const uint y)
    {
        return *m_siteReactions(x, y);
    }

    const DissolutionDeposition &surfaceReaction(const uint x, const uint y) const
    {
        return *m_siteReactions(x, y);
    }

    void setMu(const double gamma);

    const set_type &changedSurfaceSites() const
    {
        return m_changedSurfaceSites;
    }

    void updateConcentrationBoundaryIfOnBoundary(const uint x, const uint y);

    const vector<ConcentrationBoundaryReaction*> &concentrationBoundaryReactions() const
    {
        return m_concentrationBoundaryReactions;
    }

private:

    bool m_heights_set;

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

    vector<ConcentrationBoundaryReaction*> m_concentrationBoundaryReactions;

    std::vector<DissolutionDeposition*> m_affectedSurfaceReactions;

    set_type m_changedSurfaceSites;

    // KMCSolver interface
private:
    void initializeSolver();


    // Event interface
public:
    void execute();
    void initialize();
};
