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
    }

    const uint &length() const
    {
        return m_length;
    }

    const double &alpha() const
    {
        return m_alpha;
    }

    const ivec &heights() const
    {
        return m_heights;
    }

    const int &height(const uint site) const
    {
        return m_heights(site);
    }

    uint nNeighbors(const uint site) const;

    uint leftSite(const uint site, const uint n = 1) const;

    uint rightSite(const uint site, const uint n = 1) const;

    Reaction &reaction(const uint site) const
    {
        return *m_siteReactions(site);
    }

private:

    const uint m_length;

    const double m_alpha;

    ivec m_heights;

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
