#pragma once

#include "sosreaction.h"

using namespace kMC;

class SOSDiffusionReaction : public SOSReaction
{
public:
    SOSDiffusionReaction(SOSSolver &solver, const uint x0, const uint y0, const int z0);
    ~SOSDiffusionReaction();

    bool isSurfaceSite(const uint x, const uint y, const int z) const;

    const int &z() const
    {
        return m_z;
    }

    void getRandomDiffusionPath(int &dx, int &dy, int &dz) const;

private:

    int m_z;

    // Reaction interface
public:
    bool isAllowed() const;
    void executeAndUpdate();
    double rateExpression();
};

