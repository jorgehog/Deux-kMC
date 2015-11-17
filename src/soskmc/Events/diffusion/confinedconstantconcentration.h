#pragma once

#include "constantconcentration.h"


class ConfinedConstantConcentration : public ConstantConcentration
{
public:
    ConfinedConstantConcentration(SOSSolver &solver);

private:
    double m_V0;
    double m_N0;

    double m_deltaSum;

    // kMC::Observer interface
public:
    void notifyObserver(const Subjects &subject);
    void initializeObserver(const Subjects &subject);

    // Event interface
public:
    void reset();
};

