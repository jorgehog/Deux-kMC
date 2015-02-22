#pragma once

#include <kMC>
#include <armadillo>

using namespace kMC;
using namespace arma;

class SolidOnSolidSolver : public KMCSolver
{
public:
    SolidOnSolidSolver(const uint length, const double alpha);
    ~SolidOnSolidSolver();

    void registerHeightChange(const uint x, const int value)
    {
        m_heights(x) += value;
        m_heights.save("/tmp/heighmap.arma");
    }

    const uint &length() const
    {
        return m_length;
    }

    const double &alpha() const
    {
        return m_alpha;
    }

    uint nNeighbors(const uint site) const;

    uint leftSite(const uint site, const uint n = 1) const;

    uint rightSite(const uint site, const uint n = 1) const;

    bool connectedLeft(const uint leftSite, const int myHeight) const
    {
        return m_heights(leftSite) >= myHeight;
    }

    bool connectedRight(const uint rightSite, const int myHeight) const
    {
        return m_heights(rightSite) >= myHeight;
    }

    Reaction &reaction(const uint site) const
    {
        return *m_siteReactions(site);
    }

private:

    const uint m_length;

    const double m_alpha;

    vec m_heights;

    field<Reaction*> m_siteReactions;

    // KMCSolver interface
public:
    uint numberOfReactions() const;

private:
    Reaction *getReaction(const uint n) const;
};

class SolidOnSolidReaction : public Reaction
{
public:
    SolidOnSolidReaction(const uint x, SolidOnSolidSolver &system) :
        Reaction(),
        m_x(x),
        m_solver(system)
    {

    }

    SolidOnSolidSolver &system() const
    {
        return m_solver;
    }

    const uint &x() const
    {
        return m_x;
    }

    uint nNeighbors() const
    {
        return m_solver.nNeighbors(m_x);
    }

private:
    const uint m_x;

    SolidOnSolidSolver &m_solver;

};

class DiffusionDeposition : public SolidOnSolidReaction
{
public:
    using SolidOnSolidReaction::SolidOnSolidReaction;

    // Reaction interface
public:
    bool isAllowed() const;
    void executeAndUpdate();
    double rateExpression();
};














class TestReaction : public Reaction
{
public:

    TestReaction(const uint x, SolidOnSolidSolver &system) :
        Reaction(),
        m_x(x),
        m_system(system)
    {

    }

    bool isAllowed() const
    {
        return true;
    }

    void executeAndUpdate()
    {
        m_system.registerHeightChange(m_x, +1);
    }

    double rateExpression()
    {
        return 1 + sin(m_x/double(m_system.length())*2*datum::pi);
    }

private:

    const uint m_x;

    SolidOnSolidSolver &m_system;

};
