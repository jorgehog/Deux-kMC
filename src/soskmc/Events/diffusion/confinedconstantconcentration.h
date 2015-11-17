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

    // HeightObserver interface
public:
    void notifyObserver();
    void initializeObserver();

    // Event interface
public:
    void reset();
};

