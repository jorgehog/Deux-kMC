#include "pathfinder.h"

using namespace Tests;


const SearchNode *PathFinder::FindPathReversed(const int x0, const int y0, const int z0, const int x1, const int y1, const int z1)
{
    nSearchNodes = 0;
    SearchNode* startNode = makeNewSearchNode(x0, y0, z0, 0, 0, nullptr);

    MinHeap openList;
    openList.Add(startNode);

    int sx = m_world.Right();
    int sy = m_world.Top();
    m_brWorld[x0 + (y1 + z1 * sy) * sx] = true;

    while (openList.HasNext())
    {
        SearchNode* current = openList.extractFirst();

        const int r2 = current->GetDistanceSquared(x1, y1, z1);
        if (r2 <= 3)
        {
            SearchNode *result = makeNewSearchNode(x1, y1, z1, current->pathCost + r2, current->cost + r2, current);
            return result;
        }

        for (int dx = -1; dx <= 1; ++dx)
        {
            for (int dy = -1; dy <= 1; ++dy)
            {
                for (int dz = -1; dz <= 1; ++dz)
                {
                    if (dx == dy && dy == dz && dz == 0)
                    {
                        continue;
                    }

                    const int newX = current->x + dx;
                    const int newY = current->y + dy;
                    const int newZ = current->z + dz;

                    const int brWorldIdx = newX + (newY + newZ * sy) * sx;

                    if (brWorldIdx >= signed(m_brWorld.size()) || brWorldIdx < 0)
                    {
                        continue;
                    }

                    if (m_world.positionIsFree(newX, newY, newZ) && m_brWorld[brWorldIdx] == false)
                    {
                        m_brWorld[brWorldIdx] = true;
                        int pathCost = current->pathCost + dx*dx + dy*dy + dz*dz;
                        int cost = pathCost + (x1 - newX)*(x1 - newX) + (y1 - newY)*(y1 - newY) + (z1 - newZ)*(z1 - newZ);
                        SearchNode* node = makeNewSearchNode(newX, newY, newZ, cost, pathCost, current);
                        openList.Add(node);
                    }

                }
            }

        }
    }

    return nullptr; //no path found
}


