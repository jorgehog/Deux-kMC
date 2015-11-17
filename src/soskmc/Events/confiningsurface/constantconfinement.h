#pragma once

#include "fixedsurface.h"

class ConstantConfinement : public FixedSurface
{
public:
    ConstantConfinement(SOSSolver &solver, const double height);
    ~ConstantConfinement();

private:

    double m_h;

    // HeightObserver interface
public:
    void notifyObserver();

    // Event interface
public:
    void initialize();
};

