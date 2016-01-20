#pragma once

#include "../../localpotential.h"

class RDLSurface;

class RDLPotential : public LocalPotential
{
public:
    RDLPotential(SOSSolver &solver, RDLSurface &surface);

private:

    const RDLSurface &m_surface;

    // LocalPotential interface
public:
    double potential(const uint x, const uint y) const;
};
