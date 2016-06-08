#pragma once

#include "constantconcentration.h"


class ConfinedConstantConcentration : public ConstantConcentration
{
public:
    ConfinedConstantConcentration(SOSSolver &solver);

    void fixConcentration(const bool value = true)
    {
        m_fix = value;
    }

private:
    double m_N0;

    double m_deltaSum;
    double m_vCorr;

    bool m_fix;

    // kMC::Observer interface
public:
    void notifyObserver(const Subjects &subject);
    void initializeObserver(const Subjects &subject);

    // Event interface
public:
    void reset();
};

