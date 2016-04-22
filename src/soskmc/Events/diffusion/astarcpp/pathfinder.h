#pragma once

#include "minheap.h"

#include "searchnode.h"

#include "world.h"

#include <vector>

using std::vector;

namespace Tests
{    
/// <summary>
/// Author: Roy Triesscheijn (http://www.roy-t.nl)
/// Class providing 3D pathfinding capabilities using A*.
/// Heaviliy optimized for speed therefore uses slightly more memory
/// On rare cases finds the 'almost optimal' path instead of the perfect path
/// this is because we immediately return when we find the exit instead of finishing
/// 'neighbour' loop.
/// </summary>

class PathFinder
{
public:
    PathFinder(const World &world, const int nMax, const int nPreAllocNodes = 1000)
        : m_world(world),
          m_brWorld(world.Right() * world.Top() * world.Back()),
          m_nMax(nMax)
    {
        m_allSearchNodes.resize(nPreAllocNodes);

        for (int n = 0; n < nPreAllocNodes; ++n)
        {
            m_allSearchNodes[n] = new SearchNode();
        }
    }

    void clearMemory()
    {
        for (SearchNode *node : m_allSearchNodes)
        {
            delete node;
        }

        m_allSearchNodes.clear();
    }

    ~PathFinder()
    {
        clearMemory();
    }

    void reset()
    {
        for (uint i = 0; i < m_brWorld.size(); ++i)
        {
            m_brWorld[i] = false;
        }
    }

    /// <summary>
    /// Method that switfly finds the best path from start to end.
    /// </summary>
    /// <returns>The starting breadcrumb traversable via .next to the end or nullptr if there is no path</returns>
    const SearchNode* findPath(const int x0, const int y0, const int z0,
                               const int x1, const int y1, const int z1)
    {
        reset();
        //note we just flip start and end here so you don't have to.
        return FindPathReversed(x1, y1, z1, x0, y0, z0);
    }

    /// <summary>
    /// Method that switfly finds the best path from start to end. Doesn't reverse outcome
    /// </summary>
    /// <returns>The end breadcrump where each .next is a step back)</returns>
private:
    const SearchNode* FindPathReversed(const int x0, const int y0, const int z0,
                                       const int x1, const int y1, const int z1);

private:

    const World &m_world;
    vector<bool> m_brWorld;

    int m_nMax;

    vector<SearchNode*> m_allSearchNodes;
    unsigned int nSearchNodes;

    inline SearchNode *makeNewSearchNode(const int &newX, const int &newY, const int &newZ, const int &cost, const int &pathCost, SearchNode *current);

};

SearchNode *PathFinder::makeNewSearchNode(const int &newX, const int &newY, const int &newZ, const int &cost, const int &pathCost, SearchNode *current)
{
    SearchNode *nsn;

    if (nSearchNodes < m_allSearchNodes.size())
    {
        nsn = m_allSearchNodes[nSearchNodes];
        nsn->set(newX, newY, newZ, cost, pathCost, current);
    }
    else
    {
        nsn = new SearchNode(newX, newY, newZ, cost, pathCost, current);
        m_allSearchNodes.push_back(nsn);
    }

    nSearchNodes++;

    return nsn;
}

}
