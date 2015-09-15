#pragma once

#include "concentrationprofile.h"

class ConstantConcentration : public ConcentrationProfile
{
public:
    ConstantConcentration(SOSSolver &solver);
    ~ConstantConcentration();
};
