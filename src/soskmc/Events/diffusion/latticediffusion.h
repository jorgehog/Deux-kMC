#pragma once

#include "diffusion.h"

#include <unordered_map>
using std::unordered_map;

//http://stackoverflow.com/questions/22572396/using-multiple-keys-with-map-unordered-map-multidimensional
struct indices
{
    indices(const int x, const int y, const int z) :
        m_x(x),
        m_y(y),
        m_z(z)
    {

    }

    indices(const SOSDiffusionReaction *r);

    const int m_x;
    const int m_y;
    const int m_z;

    bool operator==(indices const& other) const
    {
        return std::tie(m_x, m_y, m_z) ==
               std::tie(other.m_x, other.m_y, other.m_z);
    }
};

//http://stackoverflow.com/questions/8157937/how-to-specialize-stdhashkeyoperator-for-user-defined-type-in-unordered
namespace std {
  template <> struct hash<indices>
  {
    size_t operator()(const indices & i) const
    {
        std::size_t seed = 0;
        boost::hash_combine(seed, i.m_x);
        boost::hash_combine(seed, i.m_y);
        boost::hash_combine(seed, i.m_z);

        return seed;
    }
  };
}

class LatticeDiffusion : public virtual Diffusion
{
public:

    LatticeDiffusion(SOSSolver &solver);
    ~LatticeDiffusion();

    SOSDiffusionReaction *addDiffusionReactant(const uint x, const uint y, const int z, bool setRate = true);

    void removeDiffusionReactant(SOSDiffusionReaction *reaction, bool _delete = true);

    void removeDiffusionReactant(const int x, const int y, const int z, bool _delete = true);

    SOSDiffusionReaction *diffusionReaction(const int x, const int y, const int z) const;

    void clearDiffusionReactions();

    uint numberOfDiffusionReactions() const
    {
        return m_diffusionReactionsMap.size();
    }

    void deleteQueuedReactions();

    void attachToSurface(const uint x, const uint y, const int z, SOSDiffusionReaction *reaction);

    void registerAffectedAround(const uint x, const uint y, const int z);

    void registerAffectedAroundSingle(const int neighbor, const uint xi, const uint dim, const int z);

    void dumpDiffusingParticles(const uint frameNumber, const string path = "/tmp") const;

    void moveReaction(SOSDiffusionReaction *reaction, const uint x, const uint y, const int z);

    vector<SOSDiffusionReaction *> particlesSurrounding(const uint x, const uint y, const int z) const;

private:

    vector<SOSDiffusionReaction*> dep_m_diffusionReactions;

    unordered_map<indices, SOSDiffusionReaction*> m_diffusionReactionsMap;

    vector<SOSDiffusionReaction*> m_deleteQueue;


    // Event interface
public:
    void execute();

    // Diffusion interface
public:
    virtual void dump(const uint frameNumber, const string path = "/tmp") const;
    double depositionRate(const uint x, const uint y) const;
    uint dissolutionPaths(const uint x, const uint y) const;
    void executeDiffusionReaction(SOSDiffusionReaction *reaction, const int x, const int y, const int z);
    void executeConcentrationBoundaryReaction(ConcentrationBoundaryReaction *reaction);
    bool isBlockedPosition(const uint x, const uint y, const int z) const;
    double concentration() const;


    // HeightConnecter interface
public:
    void setupInitialConditions();
    void registerHeightChange(const uint x,
                              const uint y,
                              const int value,
                              std::vector<DissolutionDeposition *> &affectedSurfaceReactions,
                              const uint n);};

