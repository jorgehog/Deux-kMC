#pragma once

#include "sosreaction.h"

class concentrationBoundaryReaction : public SOSReaction
{
public:
    concentrationBoundaryReaction(const uint dim, const uint orientation, SOSSolver &solver);
    ~concentrationBoundaryReaction();

    double freeBoundaryArea();

private:
    const uint m_dim;
    const uint m_orientation;

    SOSSolver &m_solver;

    // Reaction interface
public:
    bool isAllowed() const;
    void executeAndUpdate();
    double rateExpression();
};

