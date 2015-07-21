#pragma once

#include "../kmcsolver/reaction.h"

#include <vector>

using namespace kMC;

class SOSSolver;

class SOSReaction : public Reaction
{
public:
    SOSReaction(const uint x, const uint y, SOSSolver &system) :
        Reaction(),
        m_x(x),
        m_y(y),
        m_solver(system)
    {

    }

    SOSSolver &solver() const
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

    SOSSolver &m_solver;

};

class DiffusionDeposition : public SOSReaction
{
public:
    using SOSReaction::SOSReaction;

    const double &depositionRate() const
    {
        return m_depositionRate;
    }

    const double &dissolutionRate() const
    {
        return m_diffusionRate;
    }

    void setDiffusionRate(const double newRate)
    {
        m_diffusionRate = newRate;
        changeRate(m_diffusionRate + m_depositionRate);
    }

    void setDepositionRate(const double newRate)
    {
        m_depositionRate = newRate;
        changeRate(m_diffusionRate + m_depositionRate);
    }

    double calculateDissolutionRate() const;

    double calculateDepositionRate() const;


    // Reaction interface
public:
    bool isAllowed() const;
    void executeAndUpdate();
    double rateExpression();

private:
    double m_depositionRate;
    double m_diffusionRate;
};

