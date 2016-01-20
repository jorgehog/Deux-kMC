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

protected:

    mat m_potentialValues;

    // LocalPotential interface
public:
    double potential(const uint x, const uint y) const;

    // Observer interface
public:
    void initializeObserver(const Subjects &subject);
};

