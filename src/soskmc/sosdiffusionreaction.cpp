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
    //We are always able to diffuse down since a diffusing particle never is on the surface.
    //Transitions to the surface will transform the particle to the SOS surface and remove it
    //from the diffusion instance.
    BADAssBool(!solver().isSurfaceSite(x(), y(), z()));

    bool connectedLeft, connectedRight,
            connectedTop, connectedBottom;

    solver().findConnections(x(), y(), z(), connectedLeft, connectedRight,
                    connectedBottom, connectedTop, false);

    bool connectedAbove = solver().isBlockedPosition(x(), y(), z()+1) || solver().diffusionEvent().isBlockedPosition(x(), y(), z()+1);
    bool connectedBelow = solver().diffusionEvent().isBlockedPosition(x(), y(), z()-1);

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

    if (connectedAbove)
    {
        n--;
    }

    if (connectedBelow)
    {
        n--;
    }

    return n;

}

void SOSDiffusionReaction::getDiffusionPath(const uint path, int &dx, int &dy, int &dz)
{
    bool connectedLeft, connectedRight, connectedBottom, connectedTop;

    solver().findConnections(x(), y(), z(),
                             connectedLeft,
                             connectedRight,
                             connectedBottom,
                             connectedTop,
                             false);

    BADAssBool(!solver().isSurfaceSite(x(), y(), z()));

//    //downwards is always possible (not true anymore)
//    if (path == 0)
//    {
//        dx = 0;
//        dy = 0;
//        dz = -1;
//        return;
//    }

    uint n = 0;

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

    bool connectedBelow = solver().diffusionEvent().isBlockedPosition(x(), y(), z()-1);

    if (!connectedBelow)
    {
        if (n == path)
        {
            dx = 0;
            dy = 0;
            dz = -1;
            return;
        }
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
}

bool SOSDiffusionReaction::isAllowed() const
{
    return numberOfFreePaths() != 0;
}

void SOSDiffusionReaction::executeAndUpdate()
{
    int dx, dy, dz;

    const uint nPaths = numberOfFreePaths();

    BADAss(nPaths, !=, 0);

    uint path = rng.uniform()*nPaths;

    getDiffusionPath(path, dx, dy, dz);

    executeReaction(dx, dy, dz);
}

double SOSDiffusionReaction::rateExpression()
{
    return 6.;
}

