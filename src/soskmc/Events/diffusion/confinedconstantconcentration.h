#pragma once

#include "constantconcentration.h"


class ConfinedConstantConcentration : public ConstantConcentration
{
public:
    ConfinedConstantConcentration(SOSSolver &solver);

private:
    double m_V0;
    double m_c0;

    double m_deltaSum;

    // HeightConnecter interface
public:
    void registerHeightChange(const uint x, const uint y, const int value, std::vector<DissolutionDeposition *> &affectedSurfaceReactions, const uint nAffectedSurfaceReactions);
    void setupInitialConditions();

    // Event interface
public:
    void reset();
};

