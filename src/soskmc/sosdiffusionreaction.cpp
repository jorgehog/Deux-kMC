#include "sosdiffusionreaction.h"

#include "sossolver.h"

#include "../kmcsolver/boundary/boundary.h"
#include "../soskmc/Events/diffusion/diffusion.h"

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

uint SOSDiffusionReaction::numberOfFreePaths() const
{
    bool connectedLeft, connectedRight, connectedBottom, connectedTop;

    solver().findConnections(x(), y(), z(),
                             connectedLeft,
                             connectedRight,
                             connectedBottom,
                             connectedTop);

    //We are always able to diffuse down since a diffusing particle never is on the surface.
    //Transitions to the surface will transform the particle to the SOS surface and remove it
    //from the diffusion instance.
    BADAssBool(!solver().isSurfaceSite(x(), y(), z()));

    uint n = 6;
    if (connectedLeft)
    {
        n--;
    }

    if (connectedRight)
    {
        n--;
    }

    if (connectedBottom)
    {
        n--;
    }

    if (connectedTop)
    {
        n--;
    }

    if (solver().isBlockedPosition(x(), y(), z() + 1))
    {
        n--;
    }

    return n;

}

void SOSDiffusionReaction::getRandomDiffusionPath(int &dx, int &dy, int &dz)
{
    const uint nPaths = numberOfFreePaths();

    if (nPaths == 0)
    {
        dx = 0;
        dy = 0;
        dz = 0;

        return;
    }

    bool connectedLeft, connectedRight, connectedBottom, connectedTop;

    solver().findConnections(x(), y(), z(),
                             connectedLeft,
                             connectedRight,
                             connectedBottom,
                             connectedTop);

    uint path = rng.uniform()*nPaths;

    BADAssBool(!solver().isSurfaceSite(x(), y(), z()));

    //downwards is always possible
    if (path == 0)
    {
        dx = 0;
        dy = 0;
        dz = -1;
        return;
    }

    uint n = 1;

    if (!connectedLeft)
    {
        if (n == path)
        {
            dx = -1;
            dy = 0;
            dz = 0;

            return;
        }

        n++;
    }

    if (!connectedRight)
    {
        if (n == path)
        {
            dx = 1;
            dy = 0;
            dz = 0;

            return;
        }

        n++;
    }

    if (!connectedBottom)
    {
        if (n == path)
        {
            dx = 0;
            dy = -1;
            dz = 0;

            return;
        }

        n++;
    }

    if (!connectedTop)
    {
        if (n == path)
        {
            dx = 0;
            dy = 1;
            dz = 0;

            return;
        }

        n++;
    }

    //this leaves up the only option

    dx = 0;
    dy = 0;
    dz = 1;

    BADAssBool(!solver().isBlockedPosition(x(), y(), z() + 1));

}

void SOSDiffusionReaction::executeReaction(const int dx, const int dy, const int dz)
{
    int zNew = z() + dz;

    //Derp boundaries blocked, particle present (bunching)
    uint xNew = solver().boundary(0)->transformCoordinate(x() + dx);
    uint yNew = solver().boundary(1)->transformCoordinate(y() + dy);

    solver().diffusionEvent().executeDiffusionReaction(this, xNew, yNew, zNew);
}

void SOSDiffusionReaction::removeFromSimulation()
{
    solver().removeReaction(this);
    //derp remove from offlatticelist..
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

