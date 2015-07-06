#ifndef DIFFUSION_H
#define DIFFUSION_H

#include "../solidonsolidevent.h"

class Diffusion : public SolidOnSolidEvent
{
public:
    Diffusion(SolidOnSolidSolver &solver,
              string type,
              string unit = "",
              bool hasOutput = false,
              bool storeValue = false);

    virtual ~Diffusion();

    virtual void setupInitialConditions() = 0;

    virtual double depositionRate(const uint x, const uint y) const = 0;

    virtual void registerHeightChange(const uint x, const uint y, const int delta) = 0;

};

#endif // DIFFUSION_H
