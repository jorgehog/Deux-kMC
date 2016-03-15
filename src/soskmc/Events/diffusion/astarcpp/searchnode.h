#pragma once

#include "point3d.h"

namespace Tests
{
/// <summary>
/// Author: Roy Triesscheijn (http://www.roy-t.nl)
/// Class defining BreadCrumbs used in path finding to mark our routes
/// </summary>
class SearchNode
{
public:
    int x;
    int y;
    int z;
    int cost;
    int pathCost;
    SearchNode* next = nullptr;
    SearchNode* nextListElem = nullptr;

    SearchNode();

    SearchNode(const int x, const int y, const int z,
               int cost, int pathCost, SearchNode* next);

    ~SearchNode();

    static int refCounter;

    inline void set(const int &x, const int &y, const int &z,
             const int &cost, const int &pathCost, SearchNode* next);

    string ToString() const
    {
        stringstream ss;

        ss << x << ", " << y << ", " << z;

        return ss.str();
    }

    int GetDistanceSquared(const int x1, const int y1, const int z1) const
    {
        const int dx = x - x1;
        const int dy = y - y1;
        const int dz = z - z1;

        return (dx * dx) + (dy * dy) + (dz * dz);
    }
};

void SearchNode::set(const int &x, const int &y, const int &z, const int &cost, const int &pathCost, SearchNode *next)
{
    this->x = x;
    this->y = y;
    this->z = z;
    this->cost = cost;
    this->pathCost = pathCost;
    this->next = next;
    this->nextListElem = nullptr;
}

}
