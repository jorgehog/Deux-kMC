#pragma once

#include "firstpassagecontinuum.h"


namespace Tests
{
class World;
class PathFinder;
}

struct PathFindingJazz
{
    uint xTrans;
    uint yTrans;
    int xEnd;
    int yEnd;
    int zEnd;
};


class AStarFirstPassage : public FirstPassageContinuum
{
public:
    AStarFirstPassage(SOSSolver &solver,
                      const double maxdt,
                      const int depositionBoxHalfSize,
                      const double c);
    virtual ~AStarFirstPassage();

private:
    const int m_boxSize;

    Tests::World *m_world;
    Tests::PathFinder *m_pathFinder;
    vector<PathFindingJazz*> m_pathFindingJazzes;
    uint m_nPathFinds;


    // Event interface
public:
    void execute();

    // OfflatticeMonteCarlo interface
public:
    void calculateLocalRates();
};
