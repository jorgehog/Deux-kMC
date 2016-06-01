#pragma once

#include "sossolver.h"
#include "../kmcsolver/reaction.h"

#include "observers.h"
#include "subjects.h"

using kMC::Reaction;
using kMC::Observer;
using kMC::Subjects;

class ConcentrationBoundaryReaction : public Reaction, public Observer<Subjects>
{
public:
    ConcentrationBoundaryReaction(const uint dim,
                                  const uint orientation,
                                  SOSSolver &solver,
                                  const double omega);

    ~ConcentrationBoundaryReaction();


    uint nBoundarySites() const;
    uint nBoundarySites(const int hlInt) const;

//    double freeBoundaryArea() const;
//    double topSpacing() const;
//    double topSpacing(const double hl, const int hlInt) const;

    int heightAtBoundary(const uint n) const; //!height at boundary site n

    bool pointIsOnBoundary(const uint x, const uint y) const;

    const uint &dim() const
    {
        return m_dim;
    }

    const uint &orientation() const
    {
        return m_orientation;
    }

    const uint &location() const
    {
        return m_location;
    }

    const SOSSolver &solver() const
    {
        return m_solver;
    }

    const uint &span() const
    {
        return m_span;
    }

    double _rateExpression(const uint nBoundaryParticles) const;


private:
    const uint m_dim;
    const uint m_orientation;
    const uint m_location;
    const uint m_span;

    SOSSolver &m_solver;

    const double m_c;

    vec m_boundaryHeights;
    vec m_accuBoundaryHeights;

    void insertParticle(const uint xi, const int h);

    SOSSolver &solver()
    {
        return m_solver;
    }

    // Reaction interface
public:
    bool isAllowed() const;
    void executeAndUpdate();
    double rateExpression();

    // Observer interface
public:
    void initializeObserver(const Subjects &subject);
    void notifyObserver(const Subjects &subject);
};

