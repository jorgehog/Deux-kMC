#include "sosdiffusionreaction.h"

#include "sossolver.h"

#include "../kmcsolver/boundary/boundary.h"

#include <iostream>

using std::cout;
using std::endl;

SOSDiffusionReaction::SOSDiffusionReaction(SOSSolver &solver,
                                           const uint x0,
                                           const uint y0,
                                           const int z0) :
    SOSReaction(x0, y0, solver),
    m_z(z0)
{

}

SOSDiffusionReaction::~SOSDiffusionReaction()
{

}

bool SOSDiffusionReaction::isSurfaceSite(const uint x, const uint y, const int z) const
{
    return solver().height(x, y) == z - 1;
}

void SOSDiffusionReaction::getRandomDiffusionPath(int &dx, int &dy, int &dz) const
{
    uint dim = rng.uniform()*3;

    int direction;

    if (rng.uniform() < 0.5)
    {
        direction = 1;
    }

    else
    {
        direction = -1;
    }

    switch (dim)
    {
    case 0:
        dx = 0;
        dy = 0;
        dz = direction;
        break;

    case 1:
        dx = 0;
        dy = direction;
        dz = 0;
        break;

    case 2:
        dx = direction;
        dy = 0;
        dz = 0;
        break;
    }
}



bool SOSDiffusionReaction::isAllowed() const
{
    return true;
}

void SOSDiffusionReaction::executeAndUpdate()
{
    int dx, dy, dz;

    getRandomDiffusionPath(dx, dy, dz);

    int zNew = z() + dz;

    int roof = 100000;
    if (zNew > roof)
    {
        //Derp boundary crossed... fix this
    }

    //Derp boundaries blocked
    uint xNew = solver().boundary(0)->transformCoordinate(x() + dx);
    uint yNew = solver().boundary(0)->transformCoordinate(x() + dy);

    if (isSurfaceSite(xNew, yNew, zNew))
    {

    }
}

double SOSDiffusionReaction::rateExpression()
{
    return 1/6.;
}

