#pragma once

#include "sosreaction.h"

#include <BADAss/badass.h>

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

    uint calculateNumberOfFreePaths() const;

    void setNumberOfFreePaths();

    const uint &numberOfFreePaths() const
    {
        BADAssEqual(m_numberOfFreePaths, calculateNumberOfFreePaths());
        return m_numberOfFreePaths;
    }

    void getDiffusionPath(const uint path, int &dx, int &dy, int &dz);

    void executeReaction(const int dx, const int dy, const int dz);

private:

    int m_z;

    uint m_numberOfFreePaths;

    void removeFromSimulation();

    // Reaction interface
public:
    bool isAllowed() const;
    void executeAndUpdate();
    void affectedUpdateRule();
    double rateExpression();
};

