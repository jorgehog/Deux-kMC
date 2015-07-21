#include "sosdiffusionreaction.h"

#include <iostream>

using std::cout;
using std::endl;

SOSDiffusionReaction::SOSDiffusionReaction(const uint x0, const uint y0, const int z0)
{
    cout << "derp" << endl;
}

SOSDiffusionReaction::~SOSDiffusionReaction()
{

}



bool SOSDiffusionReaction::isAllowed() const
{
    return true;
}

void SOSDiffusionReaction::executeAndUpdate()
{
    cout << "derp" << endl;
}

double SOSDiffusionReaction::rateExpression()
{
    cout << "derp" << endl;
}
