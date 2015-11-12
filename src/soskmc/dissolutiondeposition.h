#pragma once

#include "sosreaction.h"

class DissolutionDeposition : public SOSReaction
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
    bool isAllowed() const {return true;}
    void executeAndUpdate();
    double rateExpression();

private:
    double m_depositionRate;
    double m_diffusionRate;
};
