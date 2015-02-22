#pragma once

#include <kMC>

using namespace kMC;

class SolidOnSolidSolver;

class SolidOnSolidReaction : public Reaction
{
public:
    SolidOnSolidReaction(const uint x, const uint y, SolidOnSolidSolver &system) :
        Reaction(),
        m_x(x),
        m_y(y),
        m_solver(system)
    {

    }

    SolidOnSolidSolver &solver() const
    {
        return m_solver;
    }

    const uint &x() const
    {
        return m_x;
    }

    const uint &y() const
    {
        return m_y;
    }

    uint nNeighbors() const;

private:
    const uint m_x;
    const uint m_y;

    SolidOnSolidSolver &m_solver;

};

class DiffusionDeposition : public SolidOnSolidReaction
{
public:
    using SolidOnSolidReaction::SolidOnSolidReaction;

    const double &depositionRate() const
    {
        return m_depositionRate;
    }

    const double &diffusionRate() const
    {
        return m_diffusionRate;
    }

    // Reaction interface
public:
    bool isAllowed() const;
    void executeAndUpdate();
    double rateExpression();

private:
    double m_depositionRate;
    double m_diffusionRate;

};

