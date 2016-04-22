#include "searchnode.h"

using namespace Tests;

SearchNode::SearchNode() :
    n(0),
    next(nullptr),
    nextListElem(nullptr)
{

}

SearchNode::SearchNode(const int x, const int y, const int z, int cost, int pathCost, SearchNode *next) :
    x(x),
    y(y),
    z(z),
    n(next == nullptr ? 0 : next->n + 1),
    cost(cost),
    pathCost(pathCost),
    next(next)
{

}

SearchNode::~SearchNode()
{

}


