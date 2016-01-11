#pragma once

#include <cstring>
#include <MicroPather/micropather.h>

#include <armadillo>

using micropather::Graph;
using micropather::StateCost;
using micropather::MicroPather;

class PathFinder: public Graph
{
public:

    PathFinder(const int N, const int M, const int L);

    ~PathFinder();

    void markBlockedPosition(const int x,
                             const int y,
                             const int z,
                             const int isBlocked = 1)
    {
        m_blocked(x, y, z) = isBlocked;
    }

    void getXYZ(int &x, int &y, int &z, void *state) const;

    void* area(const int x, const int y, const int z)
    {
        return (void*) (intptr_t)(z*m_N*m_M + y*m_N + x);
    }

    bool blocked(const int x, const int y, const int z)
    {
        return m_blocked(x, y, z) == 1;
    }

    int solve(std::vector<void *> &path, float &cost, const int x0, const int y0, const int z0, const int x1, const int y1, const int z1);

    void reset();

    // Graph interface
public:
    float LeastCostEstimate(void *stateStart, void *stateEnd);

    void AdjacentCost(void *state, MP_VECTOR<micropather::StateCost> *adjacent);

    void PrintStateInfo(void *state);

private:

    const int m_N;
    const int m_M;
    const int m_L;

    arma::icube m_blocked;
    arma::fvec m_costs;

    MicroPather *m_pather;
};

