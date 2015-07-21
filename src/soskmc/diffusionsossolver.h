#pragma once

#include "sossolver.h"

class DiffusionSOSSolver : public SOSSolver
{
public:

    using SOSSolver::SOSSolver;

    virtual ~DiffusionSOSSolver();

private:


    // KMCSolver interface
public:
    uint numberOfReactions() const;
    Reaction *getReaction(const uint n) const;

protected:
    void initializeSolver();

};

