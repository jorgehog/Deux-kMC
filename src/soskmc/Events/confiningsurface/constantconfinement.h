#pragma once

#include "fixedsurface.h"

class ConstantConfinement : public FixedSurface
{
public:
    ConstantConfinement(SOSSolver &solver, const double height);
    ~ConstantConfinement();

private:

    double m_h;

    // kMC::Observer interface
public:
    void notifyObserver(const Subjects &subject);

    // Event interface
public:
    void initialize();
    void execute();
};

