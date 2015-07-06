#ifndef CONSTANTCONCENTRATION_H
#define CONSTANTCONCENTRATION_H

#include "diffusion.h"

class ConstantConcentration : public Diffusion
{
public:
    ConstantConcentration(SolidOnSolidSolver &solver);

    ~ConstantConcentration();


    // Event interface
public:
    void execute()
    {
        //pass
    }

    // Diffusion interface
public:
    double depositionRate(const uint x, const uint y) const
    {
        (void) x;
        (void) y;

        return 1.0;
    }

    void registerHeightChange(const uint x, const uint y, const int delta)
    {
        (void) x;
        (void) y;
        (void) delta;

        //pass
    }

    void setupInitialConditions()
    {
        //pass
    }
};

#endif // CONSTANTCONCENTRATION_H
