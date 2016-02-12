#pragma once

#include "localpotential.h"
#include "observers.h"
#include "subjects.h"

using kMC::Subject;
using kMC::Observer;
using kMC::Subjects;


#include <armadillo>

using arma::mat;


class LocalCachedPotential : public LocalPotential, public Observer<Subjects>
{
public:
    LocalCachedPotential(SOSSolver &solver);

    virtual double potentialFunction(const uint x, const uint y) const = 0;

    const double &potentialValue(const uint x, const uint y) const
    {
        return m_potentialValues(x, y);
    }

    double sum() const;

protected:

    mat m_potentialValues;

    // LocalPotential interface
public:
    double potential(const uint x, const uint y) const;

    // Observer interface
public:
    void initializeObserver(const Subjects &subject);
};

