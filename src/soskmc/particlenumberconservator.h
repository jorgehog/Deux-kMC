#pragma once

#include "observers.h"

#include "subjects.h"

using kMC::Observer;
using kMC::Subjects;

class SOSSolver;

class ParticleNumberConservator : public Observer<Subjects>
{
public:
    ParticleNumberConservator(SOSSolver &solver);

private:

    SOSSolver &m_solver;
    uint m_targetN;

    // Observer interface
public:
    void initializeObserver(const Subjects &subject);
    void notifyObserver(const Subjects &subject);
};
