#pragma once

#include "constantconcentration.h"


class ConfinedConstantConcentration : public ConstantConcentration
{
public:
    ConfinedConstantConcentration(SOSSolver &solver);

    double newConcentration() const;

private:
    double m_V0;
    double m_c0;

    double m_deltaSum;

    double m_currentVolume;


    // HeightConnecter interface
public:
    void registerHeightChange(const uint x, const uint y, const int value, std::vector<DissolutionDeposition *> &affectedSurfaceReactions, const uint nAffectedSurfaceReactions);
    void setupInitialConditions();

    // Event interface
public:
    void reset();
};

