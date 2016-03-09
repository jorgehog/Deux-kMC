#pragma once

#include "sosreaction.h"

class SurfaceReaction : public SOSReaction
{
public:
    SurfaceReaction(const uint x, const uint y, SOSSolver &solver);

    const double &depositionRate() const
    {
        return m_depositionRate;
    }

    const double &escapeRate() const
    {
        return m_escapeRate;
    }

    void setEscapeRate(const double newRate)
    {
        m_escapeRate = newRate;
        changeRate(m_escapeRate + m_depositionRate);
    }

    void setDepositionRate(const double newRate)
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

    void freeze();

    void smelt()
    {
        m_isFrozen = false;
    }

    // Reaction interface
public:
    bool isAllowed() const {return true;}
    void executeAndUpdate();
    double rateExpression();

private:
    double m_depositionRate;
    double m_escapeRate;

    bool m_isFrozen;
    int m_freezeHeight;
};
