#include "pathfinder.h"

#include <iostream>

using std::cout;
using std::endl;


PathFinder::PathFinder(const int N, const int M, const int L) :
    Graph(),
    m_N(N),
    m_M(M),
    m_L(L),
    m_blocked(N, M, L, arma::fill::zeros),
    m_costs(3),
    m_pather(new MicroPather(this, N*M*L, 26))
{
    m_costs(0) = 1;
    m_costs(1) = sqrt(2);
    m_costs(2) = sqrt(3);
}

PathFinder::~PathFinder()
{

}

void PathFinder::getXYZ(int &x, int &y, int &z, void *state) const
{
    intptr_t index = (intptr_t)state;

    z = index/(m_N*m_M);
    const int xyIndex = index - z*(m_N*m_M);

    y = xyIndex/m_N;
    x = xyIndex - y*m_N;
}

int PathFinder::solve(std::vector<void *> path,
                      float &cost,
                      const int x0, const int y0, const int z0,
                      const int x1, const int y1, const int z1)
{
    m_pather->Reset();
    return m_pather->Solve(area(x0, y0, z0), area(x1, y1, z1), &path, &cost);
}

float PathFinder::LeastCostEstimate(void *stateStart, void *stateEnd)
{
    int x0, y0, z0;
    int x1, y1, z1;

    getXYZ(x0, y0, z0, stateStart);
    getXYZ(x1, y1, z1, stateEnd);

    int dx2 = (x0 - x1)*(x0 - x1);
    int dy2 = (y0 - y1)*(y0 - y1);
    int dz2 = (z0 - z1)*(z1 - z1);

    return sqrt(dx2 + dy2 + dz2);
}

void PathFinder::PrintStateInfo(void *state)
{
    int x, y, z;

    getXYZ(x, y, z, state);

    cout << x << ", " << y << ", " << z << endl;
}

void PathFinder::AdjacentCost(void *state, MP_VECTOR<micropather::StateCost> *adjacent)
{
    int x, y, z;
    getXYZ(x, y, z, state);

    for (int dx = -1; dx <= 1; ++dx)
    {
        for (int dy = -1; dy <= 1; ++dy)
        {
            for (int dz = -1; dz <= 1; ++dz)
            {
                if (dx == dy && dy == dz && dz == 0)
                {
                    continue;
                }

                const int nx = x + dx;
                const int ny = y + dy;
                const int nz = z + dz;


                if (nx < 0 || nx >= m_N || ny < 0 || ny >= m_M || nz < 0 || nz >= m_L)
                {
                    continue;
                }

                if (m_blocked(nx, ny, nz))
                {
                    continue;
                }

                adjacent->push_back({area(nx, ny, nz), m_costs(dx*dx + dy*dy + dz*dz - 1)});
            }
        }
    }
}

