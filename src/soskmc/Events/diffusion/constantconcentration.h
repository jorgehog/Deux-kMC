#pragma once

#include "diffusion.h"

class ConstantConcentration : public Diffusion
{
public:
    ConstantConcentration(SOSSolver &solver);

    ~ConstantConcentration();

    static double constantDepositionRate()
    {
        return 1.0;
    }

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

        return constantDepositionRate();
    }

    uint dissolutionPaths(const uint x, const uint y) const
    {
        (void) x;
        (void) y;

        return 1u;
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

    void executeDiffusionReaction(SOSDiffusionReaction *reaction, const int x, const int y, const int z)
    {
        (void) reaction;
        (void) x;
        (void) y;
        (void) z;

    }

    void executeConcentrationBoundaryReaction(ConcentrationBoundaryReaction *reaction)
    {
        (void) reaction;
    }

    bool isBlockedPosition(const uint x, const uint y, const int z) const
    {
        (void) x;
        (void) y;
        (void) z;

        return false;
    }

};
