#pragma once

#include "../kmcsolver/reaction.h"

using namespace kMC;

class SOSDiffusionReaction : public Reaction
{
public:
    SOSDiffusionReaction(const uint x0, const uint y0, const int z0);
    ~SOSDiffusionReaction();

private:

    // Reaction interface
public:
    bool isAllowed() const;
    void executeAndUpdate();
    double rateExpression();
};

