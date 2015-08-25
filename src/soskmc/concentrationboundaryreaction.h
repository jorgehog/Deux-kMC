#pragma once

#include "sossolver.h"
#include "../kmcsolver/reaction.h"

class ConcentrationBoundaryReaction : public kMC::Reaction
{
public:
    ConcentrationBoundaryReaction(const uint dim, const uint orientation, SOSSolver &solver);
    ~ConcentrationBoundaryReaction();

    double freeBoundaryArea() const;

    double dh(const uint n) const;

    const int &heightAtBoundary(const uint n) const; //!height at boundary site n

    void getFreeBoundarSite(const uint n, uint &xi, int &z) const;

    double topFilling() const;

    uint freeBoundarySites() const;

    bool pointIsOnBoundary(const uint x, const uint y) const;

    const uint &dim() const
    {
        return m_dim;
    }

    const uint &orientation() const
    {
        return m_orientation;
    }

    const uint &location() const
    {
        return m_location;
    }

    const SOSSolver &solver() const
    {
        return m_solver;
    }

    const uint &span() const
    {
        return m_span;
    }

    double _rateExpression() const;


private:
    const uint m_dim;
    const uint m_orientation;
    const uint m_location;
    const uint m_span;

    SOSSolver &m_solver;

    // Reaction interface
public:
    bool isAllowed() const;
    void executeAndUpdate();
    double rateExpression();
};

