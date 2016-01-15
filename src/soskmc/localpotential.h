#pragma once

#include <sys/types.h>

class SOSSolver;

class LocalPotential
{
public:

    //this does not add the potential to the solver.
    LocalPotential(SOSSolver &solver);

    virtual double energy(const uint x, const uint y) const = 0;

protected:

    SOSSolver &solver() const
    {
        return m_solver;
    }

private:

    SOSSolver &m_solver;

};

