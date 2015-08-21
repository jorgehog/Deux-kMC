#pragma once

#include "../solidonsolidevent.h"

class SOSDiffusionReaction;

class Diffusion : public SolidOnSolidEvent
{
public:
    Diffusion(SOSSolver &solver,
              string type,
              string unit = "",
              bool hasOutput = false,
              bool storeValue = false);

    virtual ~Diffusion();

    virtual void dump(const uint frameNumber) const;

    virtual void setupInitialConditions() = 0;

    virtual double depositionRate(const uint x, const uint y) const = 0;

    virtual void registerHeightChange(const uint x, const uint y, const int delta) = 0;

    virtual void executeDiffusionReaction(SOSDiffusionReaction *reaction, const uint x, const uint y, const int z) = 0;

    virtual bool isBlockedPosition(const uint x, const uint y, const int z) const = 0;

    virtual void insertDiffusingParticle(const double x, const double y, const double z) = 0;

};
