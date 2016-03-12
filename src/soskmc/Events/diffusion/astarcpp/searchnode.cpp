#include "searchnode.h"

using namespace Tests;

int SearchNode::refCounter = 0;


SearchNode::SearchNode() :
    next(nullptr),
    nextListElem(nullptr)
{
    refCounter++;
}

SearchNode::SearchNode(const int x, const int y, const int z, int cost, int pathCost, SearchNode *next) :
    x(x),
    y(y),
    z(z),
    cost(cost),
    pathCost(pathCost),
    next(next)
{
    refCounter++;
}

SearchNode::~SearchNode()
{
    refCounter--;
}


