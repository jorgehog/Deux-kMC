#pragma once

#include "../kmcsolver/kmcsolver.h"
#include "../kmcsolver/boundary/boundary.h"
#include "observers.h"
#include "subjects.h"

#include <armadillo>

#include <unordered_set>
#include <boost/functional/hash.hpp>

using namespace kMC;
using namespace arma;
using std::pair;
using std::unordered_set;

class SurfaceReaction;
class ConcentrationBoundaryReaction;
class ConfiningSurface;
class Diffusion;

enum class ChangeTypes
{
    Single,
    Double
};

struct CurrentSurfaceChange
{
    uint x;
    uint y;
    uint x1;
    uint y1;
    ChangeTypes type;
    int value;
};

class SOSSolver : public KMCSolver, public Subject<Subjects>, public Observer<Subjects>
{
    typedef pair<uint, uint> pair_type;
    typedef unordered_set<pair_type, boost::hash<pair_type>> set_type;

public:
    SOSSolver(const uint length,
              const uint width,
              const double alpha,
              const double gamma);

    virtual ~SOSSolver();

    void registerHeightChange(const uint x, const uint y, const int value);

    void registerSurfaceTransition(const uint x0, const uint y0, const int x1, const int y1);

    void registerChangedSite(const uint x, const uint y);

    void registerChangedAround(const uint x, const uint y);

    void setNNeighbors(const uint x, const uint y);

    void setHeight(const uint x, const uint y, const int value, const bool iteratively = true);

    void setHeights(const imat &newheights, const bool iteratively = true);

    const uint &length() const
    {
        return m_length;
    }

    const uint &width() const
    {
        return m_width;
    }

    uint minimumDimension() const
    {
        return length() < width() ? length() : width();
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
        return m_surfaceDim;
    }

    const double &alpha() const
    {
        return m_alpha;
    }

    const double &c0() const
    {
        return m_c0;
    }

    const double &gamma() const
    {
        return m_gamma;
    }

    double concentration() const
    {
        return m_concentration;
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

    bool confiningSurfaceIsSet() const
    {
        return m_confiningSurfaceEvent != nullptr;
    }

    ConfiningSurface &confiningSurfaceEvent() const
    {
        return *m_confiningSurfaceEvent;
    }

    Diffusion &diffusionEvent() const
    {
        return *m_diffusionEvent;
    }

    double volume() const;

    double freeVolume() const
    {
        return volume() - area();
    }

    double confinementEnergy(const uint x, const uint y) const;

    uint calculateNNeighbors(const uint x, const uint y, const int h) const;

    uint calculateNNeighbors(const uint x, const uint y) const
    {
        return calculateNNeighbors(x, y, height(x, y));
    }

    uint numberOfSurroundingSites(const uint x, const uint y);

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

    void getRandomSolutionSite(const uint x, const uint y, const int height,
                               int &dx, int &dy, int &dz) const;

    void getRandomSolutionSite(const uint x, const uint y,
                               int &dx, int &dy, int &dz) const
    {
        getRandomSolutionSite(x, y, height(x, y), dx, dy, dz);
    }

    int topSite(const uint x, const uint y, const int z, const uint n = 1) const;

    int bottomSite(const uint x, const uint y, const int z, const uint n = 1) const;

    int leftSite(const uint x, const uint y, const int z, const uint n = 1) const;

    int rightSite(const uint x, const uint y, const int z, const uint n = 1) const;


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

    void boundaryLatticeTransform(int &xTrans, int &yTrans, const int x, const int y, const int z) const;
    void boundaryContinousTransform(double &xTrans, double &yTrans, const double x, const double y, const double z) const;

    int boundaryLatticeTransformSingle(const int x, const int y, const int z, uint dim, const int shift = 0) const;
    double boundaryContinousTransformSingle(const double x, const double y, const double z, uint dim, const double shift = 0) const;

    void addConcentrationBoundary(const uint dim, const Boundary::orientations orientation);

    bool isBlockedPosition(const double x, const double y, const double z) const;

    bool isOutsideBoxSingle(const int x, const uint dim) const;

    bool isOutsideBox(const int x, const int y) const;

    bool isSurfaceSite(const uint x, const uint y, const int z) const;

    SurfaceReaction &surfaceReaction(const uint x, const uint y)
    {
        return *m_surfaceReactions(x, y);
    }

    const SurfaceReaction &surfaceReaction(const uint x, const uint y) const
    {
        return *m_surfaceReactions(x, y);
    }

    void setGamma(const double gamma);

    void setConcentration(const double concentration);

    void shiftConcentration(const double dc);

    const set_type &changedSurfaceSites() const
    {
        return m_changedSurfaceSites;
    }

    const vector<ConcentrationBoundaryReaction*> &concentrationBoundaryReactions() const
    {
        return m_concentrationBoundaryReactions;
    }

    double closestSquareDistance(const uint x, const uint y, const int z,
                                 const double xp, const double yp, const double zp) const;

    double absSquareDistance(const uint x, const uint y, const int z,
                             const double xp, const double yp, const double zp) const;

    const double &averageHeight() const
    {
        return m_averageHeight;
    }

    double calculateVolumeCorrection() const;

    void setZeroConcentration();

    bool concentrationIsZero() const
    {
        return m_concentrationIsZero;
    }

    const CurrentSurfaceChange &currentSurfaceChange() const
    {
        return m_currentSurfaceChange;
    }


private:

    CurrentSurfaceChange m_currentSurfaceChange;

    bool m_heights_set;

    const uint m_surfaceDim;
    ConfiningSurface *m_confiningSurfaceEvent;
    Diffusion *m_diffusionEvent;

    const uint m_length;
    const uint m_width;

    const double m_alpha;
    double m_c0;
    double m_gamma;
    double m_concentration;
    double m_expGamma;
    bool m_concentrationIsZero;

    imat m_heights;
    umat m_nNeighbors;

    double m_averageHeight;


    field<SurfaceReaction*> m_surfaceReactions;

    vector<ConcentrationBoundaryReaction*> m_concentrationBoundaryReactions;

    set_type m_changedSurfaceSites;

    // Event interface
public:
    void execute();
    void initialize();

    // KMCSolver interface
public:
    double timeUnit() const;


    // Observer interface
public:
    void initializeObserver(const Subjects &subject);
    void notifyObserver(const Subjects &subject);
};

