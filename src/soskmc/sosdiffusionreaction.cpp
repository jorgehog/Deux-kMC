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

void SOSDiffusionReaction::getRandomDiffusionPath(int &dx, int &dy, int &dz)
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

void SOSDiffusionReaction::executeReaction(const int dx, const int dy, const int dz)
{
    int zNew = z() + dz;

    int roof = 100000;
    if (zNew > roof)
    {
        //Derp boundary crossed... fix this
        removeFromSimulation();
        return;
    }

    //Derp boundaries blocked, particle present (bunching)
    uint xNew = solver().boundary(0)->transformCoordinate(x() + dx);
    uint yNew = solver().boundary(1)->transformCoordinate(y() + dy);

    if (solver().isSurfaceSite(xNew, yNew, zNew))
    {
        solver().registerHeightChange(xNew, yNew, 1);
        removeFromSimulation();
    }

    else
    {
        setX(xNew);
        setY(yNew);
    }
}

void SOSDiffusionReaction::removeFromSimulation()
{
    solver().removeReaction(this);
}

bool SOSDiffusionReaction::isAllowed() const
{
    return true;
}

void SOSDiffusionReaction::executeAndUpdate()
{
    int dx, dy, dz;

    getRandomDiffusionPath(dx, dy, dz);

    executeReaction(dx, dy, dz);
}

double SOSDiffusionReaction::rateExpression()
{
    return 1/6.;
}

