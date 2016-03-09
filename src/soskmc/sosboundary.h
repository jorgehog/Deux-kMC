#pragma once

#include "../kmcsolver/boundary/boundary.h"
#include "sossolver.h"

using namespace kMC;

class SOSBoundary : public Boundary
{
public:
    SOSBoundary(SOSSolver &solver,
                const Boundary::orientations orientation) :
        Boundary(orientation),
        m_solver(solver)
    {

    }

    virtual ~SOSBoundary() {}

    const SOSSolver &solver() const
    {
        return m_solver;
    }

protected:

    SOSSolver &solver()
    {
        return m_solver;
    }

private:

    SOSSolver &m_solver;

};
