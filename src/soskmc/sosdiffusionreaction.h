#pragma once

#include "sosreaction.h"

using namespace kMC;

class SOSDiffusionReaction : public SOSReaction
{
public:
    SOSDiffusionReaction(SOSSolver &solver, const uint x0, const uint y0, const int z0);
    ~SOSDiffusionReaction();

    const int &z() const
    {
        return m_z;
    }

    void setZ(const int z)
    {
        m_z = z;
    }

    uint numberOfFreePaths() const;

    void getRandomDiffusionPath(int &dx, int &dy, int &dz);

    void executeReaction(const int dx, const int dy, const int dz);

private:

    int m_z;

    void removeFromSimulation();

    // Reaction interface
public:
    bool isAllowed() const;
    void executeAndUpdate();
    double rateExpression();
};

