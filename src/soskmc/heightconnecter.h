#pragma once

#include <vector>
#include <sys/types.h>

class DissolutionDeposition;

//!Makes a connection between a height change and the inherited class
class HeightConnecter
{
public:
    HeightConnecter();
    virtual ~HeightConnecter();

    virtual void registerHeightChange(const uint x,
                                      const uint y,
                                      const int value,
                                      std::vector<DissolutionDeposition *> &affectedSurfaceReactions,
                                      const uint nAffectedSurfaceReactions) = 0;

    virtual void setupInitialConditions() = 0;
};
