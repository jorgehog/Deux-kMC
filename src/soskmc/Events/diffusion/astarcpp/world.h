#pragma once

#include "point3d.h"

namespace Tests
{
/// <summary>
/// Author: Roy Triesscheijn (http://www.roy-t.nl)
/// Sample World class that only provides 'is free or not' information on a node
/// </summary>
class World
{
    //Note: we use Y as height and Z as depth here!
public:
    int Left() const { return 0;  }
    int Right() const { return m_sx - 2; }
    int Bottom() const {  return 0; }
    int Top() const { return m_sy - 2; }
    int Front() const { return 0; }
    int Back() const { return m_sz - 2; }

    ~World()
    {
        delete []m_worldBlocked;
    }

    /// <summary>
    /// Creates a 3D world
    /// </summary>
    World(int width, int height, int depth = 1) :
        m_sx(width + 2),
        m_sy(height + 2),
        m_sz(depth + 2),
        m_offsetIdx((0 + 1) + ((0 + 1) + (0 + 1) * m_sy) * m_sx),
        m_worldBlocked(new bool[m_sx*m_sy*m_sz])
    {
        // added 2 to compensate for the solid border around the world

        ResetBlocks();
    }

    /// <summary>
    /// Mark positions in the world als blocked (true) or unblocked (false)
    /// </summary>
    /// <param name="value">use true if you wan't to block the value</param>
public:
    void MarkPosition(const Point3D& position, bool value)
    {
        MarkPosition(position.X, position.Y, position.Z, value);
    }

    void MarkPosition(const int X, const int Y, const int Z, bool value)
    {
        m_worldBlocked[m_offsetIdx + X + (Y + Z * m_sy) * m_sx] = value;
    }

    void ResetBlocks()
    {
        for (int i = 0; i < m_sx*m_sy*m_sz; ++i)
        {
            m_worldBlocked[i] = false;
        }

        // add solid border
        for (int x = 0; x < m_sx; ++x)
            for (int y = 0; y < m_sy; ++y)
            {
                markPositionEx(x, y, 0, true);
                markPositionEx(x, y, m_sz - 1, true);
            }

        for (int y = 0; y < m_sy; ++y)
            for (int z = 0; z < m_sz; ++z)
            {
                markPositionEx(0, y, z, true);
                markPositionEx(m_sx - 1, y, z, true);
            }

        for (int z = 0; z < m_sz; ++z)
            for (int x = 0; x < m_sx; ++x)
            {
                markPositionEx(x, 0, z, true);
                markPositionEx(x, m_sy - 1, z, true);
            }
    }

private:
    int m_sx;
    int m_sy;
    int m_sz;
    int m_offsetIdx;
    bool* m_worldBlocked; //extremely simple world where each node can be free or blocked: true=blocked


    void markPositionEx(const Point3D& position, bool value)
    {
        markPositionEx(position.X, position.Y, position.Z, value);
    }

    void markPositionEx(const int X, const int Y, const int Z, bool value)
    {
        m_worldBlocked[X + (Y + Z * m_sy) * m_sx] = value;
    }

    /// <summary>
    /// Checks if a position is free or marked (and legal)
    /// </summary>
    /// <returns>true if the position is free</returns>
public:
    bool positionIsFree(const Point3D& position) const
    {
        return positionIsFree(position.X, position.Y, position.Z);
    }

    bool positionIsFree(const int X, const int Y, const int Z) const
    {
        return !m_worldBlocked[m_offsetIdx + X + (Y + Z * m_sy) * m_sx];
    }
};

}
