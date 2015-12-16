#pragma once

#include "sosreaction.h"

class SurfaceReaction : public SOSReaction
{
public:
    using SOSReaction::SOSReaction;

    const double &depositionRate() const
    {
        return m_depositionRate;
    }

    const double &escapeRate() const
    {
        return m_escapeRate;
    }

    void setDissolutionRate(const double newRate)
    {
        m_escapeRate = newRate;
        changeRate(m_escapeRate + m_depositionRate);
    }

    void setEscapeRate(const double newRate)
    {
        m_depositionRate = newRate;
        changeRate(m_escapeRate + m_depositionRate);
    }

    double calculateEscapeRate() const;

    double calculateDepositionRate() const;

    void getEscapePath(const uint path,
                       int &dx,
                       int &dy,
                       int &dz) const;

    // Reaction interface
public:
    bool isAllowed() const {return true;}
    void executeAndUpdate();
    double rateExpression();

private:
    double m_depositionRate;
    double m_escapeRate;
};
