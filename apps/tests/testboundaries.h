#pragma once

#include <SOSkMC.h>
#include <gtest/gtest.h>

#include <sys/time.h>

#include "kmctestfixture.h"

#include "../../src/soskmc/sosdiffusionreaction.h"

//--gtest_filter=SOSkMCTest.boundaries
TEST_F(SOSkMCTest, boundaries)
{
    rng.initialize(1230303);

    const uint L = 3;
    const uint W = 3;
    const double alpha = 1.0;
    const double mu = 0;
    const double height = 20 + rng.uniform();

    m_solver = new SOSSolver(L, W, alpha, mu, getBoundariesFromIDs({0, 0, 1, 2}, L, W));
    m_pressureWallEvent = new FixedSurface(*m_solver, height);
    LatticeDiffusion *diffusionEvent = new LatticeDiffusion(*m_solver);
    m_diffusionEvent = diffusionEvent;
    SetUp_yo();

    primeSolver(0);

    for (uint x = 0; x < L; ++x)
    {
        for (uint y = 0; y < W; ++y)
        {
            m_solver->setHeight(x, y, 0);
        }
    }

    diffusionEvent->clearDiffusionReactions();

    vector<bool> isBlocked = {false, true, false};

    const Boundary *x0;
    const Boundary *x1;
    const Boundary *y0;
    const Boundary *y1;

    for (uint x0ID = 0; x0ID < isBlocked.size(); ++x0ID)
    {
        x0 = getBoundaryFromID(x0ID, L, Boundary::orientations::FIRST);

        for (uint x1ID = 0; x1ID < isBlocked.size(); ++x1ID)
        {
            x1 = getBoundaryFromID(x1ID, L, Boundary::orientations::LAST);

            for (uint y0ID = 0; y0ID < isBlocked.size(); ++y0ID)
            {
                y0 = getBoundaryFromID(y0ID, W, Boundary::orientations::FIRST);

                for (uint y1ID = 0; y1ID < isBlocked.size(); ++y1ID)
                {
                    y1 = getBoundaryFromID(y1ID, W, Boundary::orientations::LAST);

                    EXPECT_EQ(x0->isBlocked(-1), isBlocked.at(x0ID));
                    EXPECT_EQ(x1->isBlocked(L), isBlocked.at(x1ID));
                    EXPECT_EQ(y0->isBlocked(-1), isBlocked.at(y0ID));
                    EXPECT_EQ(y1->isBlocked(W), isBlocked.at(y1ID));

                    //towards center is never blocked
                    EXPECT_FALSE(x0->isBlocked(1));
                    EXPECT_FALSE(x1->isBlocked(1));
                    EXPECT_FALSE(y0->isBlocked(1));
                    EXPECT_FALSE(y1->isBlocked(1));

                    if (HasFailure())
                    {
                        return;
                    }

                    delete y1;
                }

                delete y0;
            }

            delete x1;
        }

        delete x0;
    }


    EXPECT_EQ(m_solver->calculateNNeighbors(0, 1), 5);
    EXPECT_EQ(m_solver->calculateNNeighbors(2, 1), 5);
    EXPECT_EQ(m_solver->calculateNNeighbors(1, 0), 5);
    EXPECT_EQ(m_solver->calculateNNeighbors(1, 2), 4);

    EXPECT_EQ(m_solver->nNeighbors(0, 1), 5);
    EXPECT_EQ(m_solver->nNeighbors(2, 1), 5);
    EXPECT_EQ(m_solver->nNeighbors(1, 0), 5);
    EXPECT_EQ(m_solver->nNeighbors(1, 2), 4);
}
