#pragma once

#include "diffusion.h"

class ConstantConcentration : public Diffusion
{
public:
    ConstantConcentration(SOSSolver &solver);

    ~ConstantConcentration();


    // Event interface
public:
    void execute()
    {
        const uint TRESH = 10000;

        if (cycle() % TRESH == 0)
        {
            Diffusion::dump(cycle()/TRESH);
        }

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

    void executeDiffusionReaction(SOSDiffusionReaction *reaction, const uint x, const uint y, const int z)
    {
        (void) reaction;
        (void) x;
        (void) y;
        (void) z;

    }

    bool isBlockedPosition(const uint x, const uint y, const int z) const
    {
        (void) x;
        (void) y;
        (void) z;

        return false;
    }
};
